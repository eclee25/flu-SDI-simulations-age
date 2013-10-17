#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/7/13
###Function:
##### 1) decrease susceptibility of adults in steps from 1 to 0 in an attempt to observe a "mild season" pattern in the incidence curves (T = 0.0643, where episize = 20%)

###Import data: urban_edges_Sarah.csv, urban_ages_Sarah.csv

###Command Line: python age_time_suscep.py
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
numsims = 800  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T = 0.0643 (total epidemic size = 20%)
T = 0.0643
# T = beta / (beta + gamma)
b = (-T * gamma)/(T - 1) # b = 0.0137
# define different adult susceptibilities
s1, s2 = 0, 1
susc_list = np.linspace(s1, s2, num=11, endpoint=True)

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
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/suscep_time_%ssims_beta%.3f_suscep%.1f-%.1f_vax0.zip' %(numsims, b, s1, s2)

###############################################
### susceptibility simulations ###
for s in susc_list:
	print "adult susceptibility for current sims:", s
	
	# create dict for susceptibilities
	# d_age_susc[str(age class code)] = susceptibility value
	d_age_susc = {}
	age_susc_list = [1, 1, 1, s, 1, 1]
	for c, susc in zip(range(1, 7), age_susc_list):
		d_age_susc[str(c)] = susc
	
	## save infection and recovery tsteps for each sim
	# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else float('nan')]
	d_save_I_tstep = defaultdict(list) 
	d_save_R_tstep = defaultdict(list) 
	
	for num in xrange(numsims):
		start = clock()
		child_rec, adult_rec, total_rec, OR_list, tot_incid_list, tot_preval_list, filt_OR_list, OR_tot, I_tstep_list, R_tstep_list = perc.episim_age_time_susc(G, d_node_age, T, gamma, d_age_susc)
		d_save_I_tstep[num] = I_tstep_list
		d_save_R_tstep[num] = R_tstep_list
		d_simresults[(s, num)] = (child_rec, adult_rec, total_rec)
		d_simOR[(s, num)] = OR_list
		d_simincid[(s, num)] = tot_incid_list
		d_simpreval[(s, num)] = tot_preval_list
		d_simOR_filt[(s, num)] = filt_OR_list
		d_simOR_tot[(s, num)] = OR_tot
		print "simtime, simnum:", clock()-start, "\t", num

	# print tsteps of infection and recovery to be able to recreate sim
	# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
	filename = 'Results/Itstep_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)
	
	filename = 'Results/Rtstep_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)

##############################################
### subset: epidemics only ###
# # subset epidemics from all results
# key = (susceptibility, simnumber), value = (child_rec, adult_rec, total_rec)
d_simepi = perc.epidemicsonly(d_simresults, size_epi)

# subset OR/incidence/prevalence values that produced epidemics
# key = (susceptibility, simnumber), value = [ORs/new cases/total cases/filtered ORs at different timesteps]
for key in d_simepi:
	d_epiOR[key] = d_simOR[key]
	d_epiincid[key] = d_simincid[key]
	d_epipreval[key] = d_simpreval[key]
	d_epiOR_filt[key] = d_simOR_filt[key]
	d_epiOR_tot[key[0]].append(d_simOR_tot[key])

# grab unique list of susceptibility values that produced at least one epidemic
susc_epi = list(set([key[0] for key in d_simepi]))

##############################################
### write dictionaries to files ###
# print epi OR values to file, one file per beta
for s in susc_epi:
	filename = 'Results/epiOR_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	pp.print_OR_time_to_file(d_epiOR, filename, s)
	pp.compress_to_ziparchive(zipname, filename)


##############################################
### pickle
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname2 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiincid_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epipreval_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/suscepi_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_filt_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_tot_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pickle.dump(d_epiOR, open(pname1, "wb"))
pickle.dump(d_epiincid, open(pname2, "wb"))
pickle.dump(d_epipreval, open(pname3, "wb"))
pickle.dump(susc_epi, open(pname4, "wb"))
pickle.dump(d_epiOR_filt, open(pname5, "wb"))
pickle.dump(d_epiOR_tot, open(pname6, "wb"))



