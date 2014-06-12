#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 6/3/14
###Function: Simulations where infected children and adults have a certain probability of staying home from school or work, respectively
#### Parameters able to be tuned: behavior applies to children or adults, number of days delay before cutting contacts, percentage of school/work contacts that are removed during (truncated) infectious period

###Import data: urban_ages_N10k_Sept2012.txt, urban_edges_N10k_Sept2012.txt, urban_place_N10k_Sept2012.txt

###Command Line: python age_time_sickbehavior.py
##############################################


### notes ###
# Ages:
# 1 = Infant, 2 = Toddler, 3 = Child, 4 = Adult, 5 = Senior, 6 = Elder (in nursing home)

# Places (edge attribute):
# F = household/family, S = school, H = hospital, M = shopping mall, W = workplace, D = daycare, E = elsehwere, P = preschool, O = nursing homes, N = neighbor

# T_critical = 0.0565868

### packages/modules ###
import csv
import zipfile
from time import clock
from collections import defaultdict
import networkx as nx


## local modules ##
import percolations as perc
import simulation_parameters as par
import pretty_print as pp

### simulation parameters ###
numsims = par.sp_numsims
size_epi = par.sp_size_epi
inf_period = par.sp_inf_period
gamma = par.sp_gamma
T = par.sp_T
b = par.sp_b
delay = par.sp_delay
cut = par.sp_cut
pop = par.sp_pop

print "Params:", numsims, size_epi, inf_period, gamma, T, b, delay, cut, pop

### data structures ###
edge_places_str = ''
d_node_age = {} # d_node_age[nodenumber] = ageclass
d_node_places = defaultdict(list) # d_node_places[nodenumber] = [place for contact 1, place for contact 2, .. in order of edgelist]

### import data ###
graph = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/urban_network_added_with_info_May24_2014/urban_edges_N10k_Sept2012.txt') # Vancouver network
graph_ages = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/urban_network_added_with_info_May24_2014/urban_ages_N10k_Sept2012.txt') # node number and age class
graph_places = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/urban_network_added_with_info_May24_2014/urban_place_N10k_Sept2012.txt') # node number and places where each contact are found

### construct age-structured network ###
ct = 1
G = nx.Graph()

for codes in graph_places:
	edge_places_str += codes.strip()
	d_node_places[str(ct)] = list(codes.strip())
	ct += 1
edge_places = list(edge_places_str)

for edge in graph:
	edge_ls = edge.strip().split(' ')
	G.add_edge(*edge_ls, place = edge_places.pop(0))

for line in graph_ages:
    new_line = line.strip().split(' ')
    node, age = new_line
    d_node_age[node] = age # node-ageclass dictionary

N = float(G.order())
print "network size:", N

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)

### ziparchive to write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.zip' %(numsims, b, delay, cut, pop)

###############################################
### sick behavior simulations ###
totaltime = clock()

## save infection and recovery tsteps for each sim
# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else 0]
d_save_I_tstep = defaultdict(list) 
d_save_R_tstep = defaultdict(list) 

for num in xrange(numsims):
	start = clock()
	total_rec, I_tstep_list, R_tstep_list = perc.episim_age_time_sickbehav(G, d_node_age, b, gamma, delay, cut, pop)
	d_save_I_tstep[num] = I_tstep_list
	d_save_R_tstep[num] = R_tstep_list
	print "simtime, simnum, episize:", clock() - start, "\t", num, "\t", total_rec

# print tsteps of infection and recovery to recreate sim
# sort order of sims so that the rows in d_save_I_tstep and d_save_R_tstep will match each other
filename = 'Results/Itstep_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.txt' %(numsims, b, delay, cut, pop)
pp.print_sorteddlist_to_file(d_save_I_tstep, filename, numsims)
pp.compress_to_ziparchive(zipname, filename)

filename = 'Results/Rtstep_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.txt' %(numsims, b, delay, cut, pop)
pp.print_sorteddlist_to_file(d_save_R_tstep, filename, numsims)
pp.compress_to_ziparchive(zipname, filename)

print "total time for sims:", clock() - totaltime
