#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/7/13
###Function:
##### explore relationship between OR and T in age-structured network

###Import data: urban_edges_Sarah.csv, urban_ages_Sarah.csv

###Command Line: python age_perc_T.py
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

###############################################
### data structures ###
d_node_age={} # d_node_age[nodenumber] = ageclass
d_simresults = {} # d_simresults[(T, simnumber)] = (number of infected children, number of infected adults)
d_simOR = {} # d_simOR[(T, simnumber)] = OR for epidemic or non-epidemic
d_simOR2 = defaultdict(list) # d_simOR2[T] = [OR sim1, OR sim2,...]
d_epiOR = defaultdict(list) # d_epiOR[T] = [OR epidemic1, OR epidemic2,...]
d_ressize = defaultdict(list) # d_ressize[T] = [simsize1, simsize2,...]

# for diagnostics
d_episize = defaultdict(list) # d_episize[T] = [episize1, episize2, etc...]
d_numepi = {} # d_numepi[cov] = number of epidemics produced from that coverage level

###############################################
### parameters ###
numsims = 50 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
Tlist = np.linspace(0, .2, num=21, endpoint=True) # probability of transmission
# Tlist = np.linspace(0.04, 0.07, num = 16, endpoint=True) # zoom in on realistic values

###############################################
### import data ###
f = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_edges_Sarah.csv') # Vancouver network
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

###############################################
### construct age-structured network ###
G=nx.Graph()
for edge in f:
    G.add_edge(*edge.strip().split(','))

for line in file2:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
		# node-ageclass dictionary
        d_node_age[node] = age 

net_size = float(G.order())
print "network size:", net_size

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)
print "number of child and adult nodes:", c_size, a_size


###############################################
### T simulations ###
for T in Tlist:
	print "T value for current simulations:", T
	# d_binlist[simnumber] = [list of 0s and 1s in node numbered index - 1 if node was infected in simnumber simulation]
	d_binlist = defaultdict(list) 
	
	# simulations
	start = clock()
	for num in xrange(numsims):
		child_rec, adult_rec, total_rec, bin_list, OR_val = perc.perc_age_gen(G, d_node_age, T)
		d_simresults[(T, num)] = (child_rec, adult_rec, total_rec)
		d_binlist[num] = bin_list
		d_simOR[(T, num)] = OR_val
	print "simtime, simnum:", clock()-start, "\t", num

	# print binary file of infecteds for each set of T simulations
	filename = '/home/elee/Documents/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/binlist_T_%ssims_T%.3f_vax0.txt' %(numsims, T)
# 	pp.print_dictlist_to_file(d_binlist, filename)

##############################################
### subset: epidemics only ###

# subset epidemics from all results
# key = (T, num), val = (child_rec, adult_rec, total_rec)
d_simepi = perc.epidemicsonly(d_simresults, size_epi) 

# subset ORs for sims that produced epidemics
# key = T, val = [OR epi1, OR epi2,...]
for key in sorted(d_simepi):
	d_epiOR[key[0]].append(d_simOR[key])

##############################################
### clean x-axis values ###

# unique T values that produced at least one epidemic
T_epi = sorted(list(set([key[0] for key in d_simepi])))

# total epidemic sizes for each T, epidemics only
for T in T_epi:
	# d_episize[T] = [episize1, episize2,...]
	d_episize[T] = [d_simepi[key][2] for key in sorted(d_simepi) if T == key[0]]

# separate total epidemic sizes and ORs for all results by T
for T in Tlist:

	# note: d_simresults[(T, simnumber)] = (child_rec, adult_rec, total_rec)
	# d_ressize[T] = [simsize1, simsize2,...]
	d_ressize[T] = [d_simresults[key][2] for key in sorted(d_simresults) if T == key[0]]
	
	# note: d_simOR[(T, simnumber)] = OR_val
	# d_simOR2[T] = [OR sim1, OR sim2,...]
	d_simOR2[T] = [d_simOR[key] for key in sorted(d_simOR) if T == key[0]]

##############################################
### write dictionaries to files ###
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_T_%ssims_T%.3f-%.3f_vax0.txt' %(numsims, min(Tlist), max(Tlist))
# pp.print_OR_to_file(d_epiOR, filename)

##############################################
### pickle
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname2 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_episize_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_ressize_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_simOR2_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/Tlist_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/Tepi_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname7 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_numepi_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname8 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_simepi_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pickle.dump(d_epiOR, open(pname1, "wb"))
pickle.dump(d_episize, open(pname2, "wb"))
pickle.dump(d_ressize, open(pname3, "wb"))
pickle.dump(d_simOR2, open(pname4, "wb"))
pickle.dump(Tlist, open(pname5, "wb"))
pickle.dump(T_epi, open(pname6, "wb"))
pickle.dump(d_numepi, open(pname7, "wb"))
pickle.dump(d_simepi, open(pname8, "wb"))

########################################################
# instead of using below code, use age_perc_T_viz.py module

# ##############################################
# ### RESULTS: OR by T w/ error bars ###
#
# # plot
# plt.errorbar(T_epi, [np.mean(d_epiOR[T]) for T in T_epi], yerr=[np.std(d_epiOR[T]) for T in T_epi], marker='o', color='black', linestyle='None')
# plt.xlabel('T')
# plt.ylabel('OR, child:adult')
# plt.xlim([0, 0.25])
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_T_%ssims_T%.3f-%.3f_vax0.png' %(numsims, min(Tlist), max(Tlist))
# plt.savefig(figname)
# plt.close()
# # plt.show()
#
# ##############################################
# ### RESULTS: OR by episize (epidemics only) ###
#
# # colormap formatting
# colors = cm.rainbow(np.linspace(0, 1, len(T_epi)))
#
# # separate information from dicts for each T
# for T, col in zip(T_epi, colors):
# 	episize = []
# 	# d_episize[T] = [episize1, episize2,...]
# 	for item in d_episize[T]:
# 		episize.append(item)
#
# 	# one plot per T
# 	# note: d_epiOR[T] = [OR epidemic1, OR epidemic2,...]
# 	lab_dummy = 'T=' + str(T)
# 	plt.scatter(episize, d_epiOR[T], marker = 'o', facecolors = 'none', edgecolors = col, s = 45, label = lab_dummy)
# plt.xlabel('epidemic size, epidemics only')
# plt.ylabel('OR, child:adult')
# plt.xlim([0, 10000])
# plt.ylim([0, 10])
# plt.legend(loc = 'upper left', prop = {'size':8})
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_size_T_%ssims_T%.3f_vax0.png' %(numsims, T)
# plt.savefig(figname)
# plt.close()
# # plt.show()
# 	
# ##############################################
# ### RESULTS: OR by episize (all results) ###
#
# # colormap formatting
# colors = cm.rainbow(np.linspace(0, 1, len(Tlist)))
#
# # separate information from dicts for each T
# for T, col in zip(Tlist, colors):
# 	# one plot per T value
# 	lab_dummy = 'T=' + str(T)
# 	plt.scatter(d_ressize[T], d_simOR2[T], marker = 'o', facecolors = 'none', edgecolors = col, s = 45, label = lab_dummy)
# plt.xlabel('epidemic size, all sims')
# plt.ylabel('OR, child:adult')
# plt.xlim([0, 10000])
# plt.ylim([0, 10])
# plt.legend(loc = 'upper center', prop = {'size':6})
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/simOR_size_T_%ssims_T%.3f_vax0.png' %(numsims, T)
# plt.savefig(figname)
# plt.close()
# # plt.show()
#
# ##############################################
# ### DIAGNOSTICS: epidemic size w/ error by T ###
#
# # plot episize by T
# plt.errorbar(T_epi, [np.mean(d_episize[T]) for T in T_epi], yerr=[np.std(d_episize[T]) for T in T_epi], marker='o', color='black', linestyle='None')
# plt.xlabel('T')
# plt.ylabel('epidemic size')
# plt.xlim([0, 0.25])
# plt.ylim([0, 10000])
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/episize_T_%ssims_T%.3f-%.3f_vax0.png' %(numsims, min(Tlist), max(Tlist))
# plt.savefig(figname)
# plt.close()
# # plt.show()
#
# ##############################################
# ### DIAGNOSTICS: number of epidemics by T ### 
# for T in T_epi:
# 	d_numepi[T] = len([d_simepi[key][2] for key in d_simepi if T == key[0]])
#
# # plot
# plt.plot(sorted(T_epi), [d_numepi[T] for T in sorted(T_epi)], marker='o', color='black')
# plt.xlabel('T')
# plt.ylabel('number of epidemics')
# plt.xlim([0, 0.25])
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/numepi_T_%ssims_T%.3f-%.3f_vax0.png' %(numsims, min(Tlist), max(Tlist))
# plt.savefig(figname)
# plt.close()
# # plt.show()



