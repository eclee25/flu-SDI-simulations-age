#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 12/23/13
###Function:
##### 1) Change child infectious periods in steps relative to infectious periods of other age groups (while infectious period in all other age groups remains at 5 days) in an attempt to observe a "mild season" pattern in the incidence curves (T = 0.0643, where episize = 20%)

### Note: Simulations run fairly slow at longer infectious periods (as expectd)

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
numsims = 800  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
inf_period = 5. # 5 days recovery for all non-children
gamma = 1/inf_period
# base values of T
# see calculation for beta within simulation script
T = 0.0643 # total epidemic size (naive, no age-dep params) = 20%
# T = 0.075 # total epidemic size (naive, no age-dep params) = 30%
# calculate base b for non-children when T = beta / (beta + gamma)
# when T = 0.0643 and gamma = 1/5, b = 0.0137
# when T = 0.075 and gamma = 1/5, b = 0.0162
b = (-T * gamma)/(T - 1) 

# define different adult infectious periods ranging from 3 to 15 days
# Cauchemez2004 reports that children have an infectious period of 3.6 days and adults have an infectious period of 3.9 days (relative relationship is 3.6/3.9 = 0.923)
# proportion of child infectious period relative to that of non-child infectious period 
r1, r2 = 3, 15
rec_list = np.linspace(r1, r2, num=9, endpoint=True) 


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
print "network size:", N

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)

### ziparchive to write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/childrecov_time_%ssims_beta%.3f_rec%.1f-%.1f_vax0.zip' %(numsims, b, r1, r2)

###############################################
### infectious period length simulations ###
totaltime = clock()

for r in rec_list:
	print "child gamma for current sims:", r
	
	# create dict for child and non-child betas, which are age-dependent on gamma, and gammas
	# d_age_rec[str(age class code)] = (beta, gamma)
	d_age_rec = {}
	age_rec_list = [gamma, gamma, 1/r, gamma, gamma, gamma] 
	# children are the third age class in d_node_age, adults are the fourth
	for ageclass, rec in zip(range(1, 7), age_rec_list):
		d_age_rec[str(ageclass)] = (b, rec)
	print "ageclass code, (beta, gamma):", d_age_rec.items()
	
	## save infection and recovery tsteps for each sim
	# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else 0]
	d_save_I_tstep = defaultdict(list) 
	d_save_R_tstep = defaultdict(list) 
	
	# timer for all sims of one recovery value
	start_all = clock()
	for num in xrange(numsims):
		start = clock()
		total_rec, I_tstep_list, R_tstep_list = perc.episim_age_time_rec(G, d_node_age, d_age_rec)
		d_save_I_tstep[num] = I_tstep_list
		d_save_R_tstep[num] = R_tstep_list
		print "simtime, simnum, episize:", clock() - start, "\t", num, "\t", total_rec
	print "simtime for %s sims for child gamma %1.1f" %(numsims, r), clock() - start_all

# print tsteps of infection and recovery to recreate sim
# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
	filename = 'Results/Itstep_childrecov_time_%ssims_beta%.3f_rec%.1f_vax0.txt' %(numsims, b, r)
	pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)
	
	filename = 'Results/Rtstep_childrecov_time_%ssims_beta%.3f_rec%.1f_vax0.txt' %(numsims, b, r)
	pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)

print "total time for sims:", clock() - totaltime

# update below code for T_avg with varying gammas
T_na = b/(b + gamma)
for r in rec_list:
	gamma_a = 1/r
	T_a = b/(b + gamma_a)
	T_avg = (T_a * a_size + T_na * (N - a_size))/N
	print "gamma_a, T_na, T_a, T_avg:", gamma_a, T_na, T_a, T_avg

