#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 11/20/13
###Function: Time-based epidemic simulations with realistic parameters for children and non-children. Non-children are considered the base age class. Realistic parameters include: individual susceptibility (sigma), individual infectivity per contact (beta), individual recovery rate (gamma = 1/infectious period). Realistic parameters come from Cauchemez, et al. (2004). A Bayesian MCMC approach to study transmission of influenza: application to household longitudinal data. Statistics in Medicine 23: 3469-3487.

###Import data: urban_edges_Sarah.csv, urban_ages_Sarah.csv

###Command Line: python age_time_realistic.py
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
# d_age_params[str(age class code)] = (sigma, beta, gamma)
d_age_params = {}
# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else 0]
d_save_I_tstep = defaultdict(list) 
d_save_R_tstep = defaultdict(list)

### simulation parameters ###
numsims = 10  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)

## realistic disease parameters (Cauchemez 2004)
# gamma = probability of recovery at each time step
# 3.6 day infectious period for children
# 3.9 day infectious period for non-children
gamma_c, gamma_nc = 1/3.6, 1/3.9  
# sigma = susceptibility of individual
sigma_c, sigma_nc = 1.15, 1.0
# infectivity = probability of infecting a contact over entire infectious period
infectivity_c, infectivity_nc = 0.48, 0.26
infectivity_cnc = infectivity_c/infectivity_nc

## T defined by epidemic size on network
# T = transmissibility of individual
# T_nc = 0.0643 # total epidemic size (naive, no age-dep params) = 20%
# T_nc = 0.075 # total epidemic size (naive, no age-dep params) = 30%
T_nc = 0.0502 # T_avg = 0.0643
T_c = T_nc * infectivity_cnc
# T = beta / (beta + gamma)
b_c = (-T_c * gamma_c)/(T_c - 1) 
b_nc = (-T_nc * gamma_nc)/(T_nc - 1)

# assign parameters to dict for ease of importing into a function
c_param_tuple = (sigma_c, b_c, gamma_c)
nc_param_tuple = (sigma_nc, b_nc, gamma_nc)
for ageclass in range(1, 7):
	d_age_params[str(ageclass)] = nc_param_tuple
d_age_params['3'] = c_param_tuple

### functions ###

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
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.zip' %(numsims, b_nc)

###############################################
### susceptibility simulations ###

totaltime = clock()
print "child: sigma, beta, gamma", c_param_tuple
print "non-child: sigma, beta, gamma", nc_param_tuple

for num in xrange(numsims):
	start = clock()
	total_rec, I_tstep_list, R_tstep_list = perc.episim_age_time_realistic(G, d_node_age, d_age_params)
	d_save_I_tstep[num] = I_tstep_list
	d_save_R_tstep[num] = R_tstep_list
	print "simtime, simnum, episize:", clock() - start, "\t", num, "\t", total_rec

# print infection and recovery tsteps to recreate sim
# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
filename = 'Results/Itstep_realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.txt' %(numsims, b_nc)
pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
pp.compress_to_ziparchive(zipname, filename)

filename = 'Results/Rtstep_realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.txt' %(numsims, b_nc)
pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
pp.compress_to_ziparchive(zipname, filename)

print "for %s sims for nonchild params:" %(numsims), nc_param_tuple, "\tsimtime", clock() - totaltime

T_avg = (T_c * sigma_c * c_size + T_nc * (N - c_size))/N
print "T_avg:", T_avg

























