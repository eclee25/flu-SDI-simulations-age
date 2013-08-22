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

## local modules ##
import percolations as perc

### data structures ###
d_node_age={} # d_node_age[nodenumber] = ageclass
d_simresults = {} # d_simresults[(beta, simnumber)] = number infected
d_episize = defaultdict(list) # epidemic sizes of epidemics only
d_simOR = defaultdict(list) # OR for each time step for all sims
d_epiOR = defaultdict(list) # OR for each time step for epidemics only

### parameters ###
numsims = 100 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
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


###############################################
### beta simulations ###
for beta in blist:
	print "beta value for current simulations:", beta
	# d_binlist[simnumber] = [list of 0s and 1s in node numbered index - 1 if node was infected in entire simnumber simulation]
	d_binlist = defaultdict(list) 
	for num in xrange(numsims):
		start = clock()
		child_rec, adult_rec, total_rec, bin_list, OR_list = perc.perc_age_time(G, d_node_age, beta, gamma)
		d_binlist[num] = bin_list
		d_simresults[(beta, num)] = (child_rec, adult_rec, total_rec)
		d_simOR[(beta, num)] = OR_list
		print "simtime, simnum:", clock()-start, "\t", num

	# print binary file of infecteds for each set of T simulations
	# order of simulations in dictionary doesn't matter
	filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/binlist_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	pp.print_dictlist_to_file(d_binlist, filename)


##############################################
### subset: epidemics only ###
# # subset epidemics from all results
# key = (beta, simnumber), value = (child_rec, adult_rec, total_rec)
d_simepi = perc.epidemicsonly(d_simresults, size_epi)

# subset OR values that produced epidemics
# key = (beta, simnumber), value = [ORs at different timesteps]
for key in d_simepi:
	d_epiOR[key] = d_simOR[key]

# grab unique list of betas that produced at least one epidemic
beta_epi = list(set([key[0] for key in d_simepi]))

##############################################
### plot OR by time for each beta value ###
# each sim is one line
for beta in beta_epi:
	pl_ls = [key for key in d_epiOR if key[0] == beta]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
	plt.plot(xrange(len(d_epiOR[key])), [1] * len(d_epiOR[key]), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, beta: ' + str(beta))
	plt.ylabel('OR, child:adult')
	plt.ylim([-3, 25])
	plt.xlim([-1, 100])
	figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_beta_time_%ssims_beta%.3f_vax0_%s.png' %(numsims, beta, key)
	plt.savefig(figname)
	plt.close()
# 	plt.show()


#
# ### plot OR by episize (epidemics only) ###
# episize, epiOR = [],[]
# for T in sorted(T_epi):
# 	for item in d_episize[T]:
# 		episize.append(item)
# 	for item in d_epiOR[T]:
# 		epiOR.append(item)
# plt.scatter(episize, epiOR, marker = 'o', facecolors='none', edgecolors='blue', s = 45, label='epidemics only')
# plt.xlabel('epidemic size')
# plt.ylabel('OR, child:adult')
# plt.legend(loc = 'upper left')
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_size_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
# plt.savefig(figname)
# # plt.show()
# plt.close()
#
# ##############################################
# ### DIAGNOSTICS: epidemic size w/ error bars ###
# # plot episize by T
# plt.errorbar(T_epi, [np.mean(d_episize[T]) for T in T_epi], yerr=[d_episize_sd[T] for T in T_epi], marker='o', color='black', linestyle='None')
# plt.xlabel('T')
# plt.ylabel('epidemic size')
# plt.xlim([-.05, 0.25])
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/episize_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
# plt.savefig(figname)
# # plt.show()
# plt.close()
#
# ### DIAGNOSTICS: number of epidemics ### 
# d_numepi ={}
# T_epi = list(set([key[0] for key in d_simepi]))
# for T in T_epi:
# 	d_numepi[T] = len([d_simepi[key][2] for key in d_simepi if T == key[0]])
#
# # plot
# plt.plot(sorted(T_epi), [d_numepi[T] for T in sorted(T_epi)], marker='o', color='black')
# plt.xlabel('T')
# plt.ylabel('number of epidemics')
# plt.xlim([-.05, 0.25])
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/numepi_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
# plt.savefig(figname)
# # plt.show()
# plt.close()
#
#
##############################################
### write dictionaries to files ###
# print epi OR values to file, one file per beta
for beta in beta_epi:
	filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	pp.print_OR_time_to_file(d_epiOR, filename, beta)



