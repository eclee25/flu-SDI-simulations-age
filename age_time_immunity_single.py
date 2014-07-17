#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 7/15/14
###Function: Conduct time-based simulations on age-structure networks where pre-existing immunity exists and is heterogeneous within the adult population. Vary single average value of immunity for subpopulation of adults with any pre-existing immunity. Choose subpopulation node IDs with a fixed random proportion of the adult population.

###Import data: 

###Command Line: python age_time_immunity_single.py
##############################################


### notes ###


### packages/modules ###
import csv
import zipfile
from time import clock
from collections import defaultdict
import networkx as nx
from random import seed

## local modules ##
import percolations as perc
import simulation_parameters as par
import pretty_print as pp

### parameters ###
seed(19)
numsims = par.sp_numsims
size_epi = par.sp_size_epi
inf_period = par.sp_inf_period
g = par.sp_gamma
T = par.sp_T
b = par.sp_b
# specific to immunity params
imm_val_ls = par.sp_immune_val_list
prop = par.sp_prop
zstring = par.sp_pstr_fixed
zstring2 = par.sp_mstr_range

### data structures ###
d_node_age = {} # d_node_age[nodenumber] = ageclass

###############################################
### import data and initialize graph ###

### import data ###
graph = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/urban_network_added_with_info_May24_2014/urban_edges_N10k_Sept2012.txt') # Vancouver network
graph_ages = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/urban_network_added_with_info_May24_2014/urban_ages_N10k_Sept2012.txt') # node number and age class

### construct age-structured network ###
ct = 1
G = nx.Graph()

for edge in graph:
	edge_ls = edge.strip().split(' ')
	G.add_edge(*edge_ls)

for line in graph_ages:
    new_line = line.strip().split(' ')
    node, age = new_line
    d_node_age[node] = age # node-ageclass dictionary

N = float(G.order())
print "network size:", N

### ziparchive to write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/immunity_time_%ssims_beta%.3f_%s_%s.zip' %(numsims, b, zstring, zstring2)

###############################################
### set pre-existing immunity conditions ###
totalsimtime = clock()
## identify nodes with pre-existing immunity (return list of adult IDs - mult methods) ## imm_nodes = list of adult IDs with any pre-existing immunity

# choose randomly, given a proportion of adults with pre-existing immunity (# tune adult proportion in equal increments) (# set adult proportion equal to the proportion of adults infected from a single simulation)
imm_nodes = perc.immune_nodes_proportion(G, d_node_age, prop)	
	
for imm_val in imm_val_ls:
	zstring3 = 'single%s' %(imm_val) # string for filename disambiguation

	## assign magnitude of pre-existing immunity to each node (return dict with node ID and pre-existing immunity value; 0 if none - mult methods) ## d_immunity_mag[node] = pre-existing immunity value

	# set a single average value of immunity
	d_immunity_mag = perc.immunity_magnitude_single(G, imm_nodes, imm_val)


	###############################################
	### pre-existing immunity simulations ###
	totaltime = clock()

	## save infection and recovery tsteps for each sim
	# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else 0]
	d_save_I_tstep = defaultdict(list) 
	d_save_R_tstep = defaultdict(list) 

	for num in xrange(numsims):
		start = clock()
		total_rec, I_tstep_list, R_tstep_list = perc.episim_age_time_imm(G, d_node_age, d_immunity_mag, b, g)
		d_save_I_tstep[num] = I_tstep_list
		d_save_R_tstep[num] = R_tstep_list
		print "simtime, simnum, episize:", clock() - start, "\t", num, "\t", total_rec

	# print tsteps of infection and recovery to recreate sim
	# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
	filename = 'Results/Itstep_immunity_time_%ssims_beta%.3f_%s_%s.txt' %(numsims, b, zstring, zstring3)
	pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)

	filename = 'Results/Rtstep_immunity_time_%ssims_beta%.3f_%s_%s.txt' %(numsims, b, zstring, zstring3)
	pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
	pp.compress_to_ziparchive(zipname, filename)

	print "total time for sims:", clock() - totaltime
	print "Params:", numsims, size_epi, inf_period, g, T, b, imm_val, prop

print "total time for sims:", clock() - totalsimtime
print "age_time_immunity_single.py complete"