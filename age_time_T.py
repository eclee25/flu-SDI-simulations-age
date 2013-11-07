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
numsims = 50  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T ranges from 0.0 to 0.2, gamma = 1/5 and T = beta / (beta + gamma)
# T1, T2 = 0.0, 0.2
# T1, T2 = 0.075, 0.075 # compare to age_time_suscep_viz
T1, T2 = 0.0643, 0.0643
b1, b2 = (-T1 * gamma)/(T1 - 1), (-T2 * gamma)/(T2 - 1) # 0, .05
blist = np.linspace(b1, b2, num=1, endpoint=True) # probability of transmission

### import data ###
graph = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_edges_Sarah.csv') # Vancouver network
graph_ages = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

### construct age-structured network ###
G=nx.Graph()
for edge in graph:
    G.add_edge(*edge.strip().split(','))

for line in graph_ages:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age # node-ageclass dictionary

net_size = float(G.order())
print "network size:", net_size

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)


### ziparchive to write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/beta_time_%ssims_beta%.3f-%.3f_vax0.zip' %(numsims, b1, b2)

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
		# 11/5/13 reduced outputs of simulation
		total_rec, I_tstep_list, R_tstep_list = perc.episim_age_time(G, d_node_age, beta, gamma)
		d_save_I_tstep[num] = I_tstep_list
		d_save_R_tstep[num] = R_tstep_list
		print "simtime, simnum, episize:", clock()-start, "\t", num, "\t", total_rec

	# print tsteps of infection and recovery to be able to recreate sim
	# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
	filename = 'Results/Itstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)
	
	filename = 'Results/Rtstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)


#
#
# ##############################################
# ### write dictionaries to files ###
# # print epi OR values to file, one file per beta
# for beta in beta_epi:
# 	filename = 'Results/epiOR_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
# 	pp.print_OR_time_to_file(d_epiOR, filename, beta)
# 	pp.compress_to_ziparchive(zipname, filename)
#
#
# ##############################################
# ### pickle
# pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
# pname2 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiincid_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
# pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epipreval_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
# pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/betaepi_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
# pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_filt_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
# pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_tot_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
# pickle.dump(d_epiOR, open(pname1, "wb"))
# pickle.dump(d_epiincid, open(pname2, "wb"))
# pickle.dump(d_epipreval, open(pname3, "wb"))
# pickle.dump(beta_epi, open(pname4, "wb"))
# pickle.dump(d_epiOR_filt, open(pname5, "wb"))
# pickle.dump(d_epiOR_tot, open(pname6, "wb"))


