#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/7/13
###Function:
##### 1) odds ratio by T while keeping track of time

###Import data: urban_edges_Sarah.csv, urban_ages_Sarah.csv

###Command Line: python age_perc.py
##############################################

### notes ###
# codebook of age class codes
# '1' - Toddlers: 0-2
# '2' - Preschool: 3-4
# '3' - Children: 5-18
# '4' - Adults: 19-64
# '5' - Seniors: 65+ (community)
# '6' - Elders: 65+ (nursing home)
# There are only 94 "elders" in the Vancouver network, and they all reside in one nursing home, so they can be combined with the seniors for analysis purposes (all_elderly).

### packages/modules ###
import networkx as nx
import random as rnd
import numpy as np
from time import clock
from collections import defaultdict
import matplotlib.pyplot as plt
import pretty_print as pp
import pickle 

## local modules ##
import percolations as perc

### data structures ###
d_node_age={} # d_node_age[nodenumber] = ageclass
d_simresults = {} # d_simresults[(beta, simnumber)] = number infected
d_episize = defaultdict(list) # epidemic sizes of epidemics only
d_simOR = defaultdict(list) # OR for each time step for all sims
d_epiOR = defaultdict(list) # OR for each time step for epidemics only
d_filt_tsteps = defaultdict(list) # start and end time steps for each simulation by which we exclude OR data on chart due to small infected numbers
d_simOR_filt = defaultdict(list) # OR for each time step for all sims where OR is nan when we want to exclude the time point due to small infected numbers
d_epiOR_filt = defaultdict(list) # OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers
d_simincid = defaultdict(list) # incidence for each time step for all sims
d_epiincid = defaultdict(list) # incidence for each time step for epidemics only
d_simpreval = defaultdict(list) # prevalence for each time step for all sims
d_epipreval = defaultdict(list) # prevalence for each time step for epidemics only
d_simOR_tot = {} # d_simOR_tot[(beta, simnumber)] = OR for entire simulation for all results
d_epiOR_tot = defaultdict(list) # d_epiOR_tot[beta] = list of ORs for all simulations that were epidemics

### parameters ###
numsims = 1000  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T ranges from 0.0 to 0.2, gamma = 1/5 and T = beta / (beta + gamma)
b1, b2 = (-0.0 * gamma)/(0 - 1), (-0.2 * gamma)/(0.2 - 1) # 0, .05
blist = np.linspace(b1, b2, num=11, endpoint=True) # probability of transmission

### import data ###
f = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_edges_Sarah.csv') # Vancouver network
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

### construct age-structured network ###
G=nx.Graph()
for edge in f:
    G.add_edge(*edge.strip().split(','))

for line in file2:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age # node-ageclass dictionary

net_size = float(G.order())
print "network size:", net_size

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)


### ziparchive to write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/beta_time_%ssims_beta%.3f-beta%.3f_vax0.zip' %(numsims, b1, b2)

###############################################
### beta simulations ###
for beta in blist:
	print "beta value for current simulations:", beta
	
	## save infection and recovery tsteps for each sim
	# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else float('nan')]
	d_save_I_tstep = defaultdict(list) 
	d_save_R_tstep = defaultdict(list) 
	
	for num in xrange(numsims):
		start = clock()
		child_rec, adult_rec, total_rec, OR_list, tot_incid_list, tot_preval_list, filt_OR_list, OR_tot, I_tstep_list, R_tstep_list = perc.episim_age_time(G, d_node_age, beta, gamma)
		d_save_I_tstep[num] = I_tstep_list
		d_save_R_tstep[num] = R_tstep_list
		d_simresults[(beta, num)] = (child_rec, adult_rec, total_rec)
		d_simOR[(beta, num)] = OR_list
		d_simincid[(beta, num)] = tot_incid_list
		d_simpreval[(beta, num)] = tot_preval_list
		d_simOR_filt[(beta, num)] = filt_OR_list
		d_simOR_tot[(beta, num)] = OR_tot
		print "simtime, simnum:", clock()-start, "\t", num

	# print tsteps of infection and recovery to be able to recreate sim
	# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
	filename = 'Results/Itstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)
	
	filename = 'Results/Rtstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)

##############################################
### subset: epidemics only ###
# # subset epidemics from all results
# key = (beta, simnumber), value = (child_rec, adult_rec, total_rec)
d_simepi = perc.epidemicsonly(d_simresults, size_epi)

# subset OR/incidence/prevalence values that produced epidemics
# key = (beta, simnumber), value = [ORs/new cases/total cases/filtered ORs at different timesteps]
for key in d_simepi:
	d_epiOR[key] = d_simOR[key]
	d_epiincid[key] = d_simincid[key]
	d_epipreval[key] = d_simpreval[key]
	d_epiOR_filt[key] = d_simOR_filt[key]
	d_epiOR_tot[key[0]].append(d_simOR_tot[key])

# grab unique list of betas that produced at least one epidemic
beta_epi = list(set([key[0] for key in d_simepi]))

##############################################
### calculate avg OR for each time point


##############################################
### write dictionaries to files ###
# print epi OR values to file, one file per beta
for beta in beta_epi:
	filename = 'Results/epiOR_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	pp.print_OR_time_to_file(d_epiOR, filename, beta)
	pp.compress_to_ziparchive(zipname, filename)


##############################################
### pickle
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname2 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiincid_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epipreval_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/betaepi_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_filt_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_tot_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pickle.dump(d_epiOR, open(pname1, "wb"))
pickle.dump(d_epiincid, open(pname2, "wb"))
pickle.dump(d_epipreval, open(pname3, "wb"))
pickle.dump(beta_epi, open(pname4, "wb"))
pickle.dump(d_epiOR_filt, open(pname5, "wb"))
pickle.dump(d_epiOR_tot, open(pname6, "wb"))


########################################################
# instead of plotting in this script, use age_perc_T_time_viz.py module

##############################################
### plot OR by time for each beta value ###
# # each sim is one line
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epiOR if key[0] == beta]
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
# 	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
# 	plt.xlabel('time step, beta: ' + str(beta))
# 	plt.ylabel('OR, child:adult')
# 	plt.ylim([-3, 15])
# 	plt.xlim([-1, 125])
# 	figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
# 	plt.savefig(figname)
# 	plt.close()
# 	plt.show()
#
# ##############################################
# ### plot incidence by time for each beta value ###
# # each sim is one line
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epiincid if key[0] == beta]
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')
# 	plt.xlabel('time step, beta: ' + str(beta))
# 	plt.ylabel('number of new cases')
# 	plt.xlim([-1, 125])
# 	figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiincid_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
# 	plt.savefig(figname)
# 	plt.close()
# 	plt.show()
#
# ##############################################
# ### plot prevalence by time for each beta value ###
# # each sim is one line
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epipreval if key[0] == beta]
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epipreval[key])), d_epipreval[key], marker = 'None', color = 'grey')
# 	plt.xlabel('time step, beta: ' + str(beta))
# 	plt.ylabel('total number of cases')
# 	plt.xlim([-1, 125])
# 	figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epipreval_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
# 	plt.savefig(figname)
# 	plt.close()
# 	plt.show()

