#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/7/13
###Function:
##### 1) increase susceptibility of children in steps from 1 to 1.5 in an attempt to observe a "mild season" pattern in the incidence curves (T = 0.0643, where episize = 20%)

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
# T_critical = 0.0565868

### packages/modules ###
import networkx as nx
import random as rnd
import numpy as np
from time import clock
from collections import defaultdict
import matplotlib.pyplot as plt
import zipfile

## local modules ##
import percolations as perc
import pretty_print as pp

### data structures ###
# d_node_age[nodenumber] = ageclass
d_node_age = {} 

### simulation parameters ###
numsims = 10  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 1/float(5) # 5 days recovery here
# T = 0.0643 # total epidemic size (naive, no age-dep params) = 20%
# T = 0.075 # total epidemic size (naive, no age-dep params) = 30%
T = 0.0620 # T_avg = 0.0643 @ sigma_c = 1.15

# T = beta / (beta + gamma)
# when T = 0.0643 and gamma = 1/5, b = 0.0137
# when T = 0.075 and gamma = 1/5, b = 0.0162
b = (-T * gamma)/(T - 1) 

# define different child susceptibilities
# Cauchemez 2004 cites child susceptibility to be 1.15 times greater than that of adults
s1, s2 = 1, 1.5
susc_list = np.linspace(s1, s2, num=6, endpoint=True)

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

N = float(G.order())
print "network size:", net_size

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)

### ziparchive to write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/suscep_time_%ssims_beta%.3f_suscep%.1f-%.1f_vax0.zip' %(numsims, b, s1, s2)

###############################################
### susceptibility simulations ###
totaltime = clock()

for s in susc_list:
	print "child susceptibility for current sims:", s
	
	# create dict for susceptibilities
	# d_age_susc[str(age class code)] = susceptibility value
	d_age_susc = {}
	# children are the third age class in d_node_age
	age_susc_list = [1, 1, s, 1, 1, 1] 
	for ageclass, susc in zip(range(1, 7), age_susc_list):
		d_age_susc[str(ageclass)] = susc
	print d_age_susc.items()
	
	## save infection and recovery tsteps for each sim
	# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else 0]
	d_save_I_tstep = defaultdict(list) 
	d_save_R_tstep = defaultdict(list) 
	
	# timer for all sims of one child susceptibility
	start_all = clock()
	for num in xrange(numsims):
		start = clock()
		total_rec, I_tstep_list, R_tstep_list = perc.episim_age_time_susc(G, d_node_age, b, gamma, d_age_susc)
		d_save_I_tstep[num] = I_tstep_list
		d_save_R_tstep[num] = R_tstep_list
		print "simtime, simnum, episize:", clock() - start, "\t", num, "\t", total_rec
	print "simtime for %s sims for child suscep %1.1f" %(numsims, s), clock() - start_all

# print tsteps of infection and recovery to recreate sim
# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
	filename = 'Results/Itstep_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)
	
	filename = 'Results/Rtstep_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)

print "total time for sims:", clock() - totaltime

# reference table: probability of infection before adjusting for susceptibility
for inf in range(52):
	print inf, 1- np.exp(-b * inf)

T_avg = (T * 1.15 * c_size + T * (N - c_size))/N
print "T_avg:", T_avg
