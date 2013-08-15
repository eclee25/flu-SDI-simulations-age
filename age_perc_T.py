#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/7/13
###Function:
##### 1) odds ratio by T
##### 2) odds ratio by vax efficacy

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
d_simresults = {} # d_simresults[(T, simnumber)] = (number of infected children, number of infected adults)
d_resincid = defaultdict(list) # d_resincid[(T, 'C' or 'A')] = [child/adult incidence simulation1, child/adult incidence simulation2,...]
d_resOR = defaultdict(list) # d_resOR[T] = [OR sim1, OR sim2,...]
d_resORsd = {} # d_resORsd[T] = sd of ORs across all simulations
d_epiincid = defaultdict(list) # d_epiincid[(T, 'C' or 'A')] = [child/adult incidence epidemic1, child/adult incidence epidemic2,...]
d_epiOR = defaultdict(list) # d_epiOR[T] = [OR epidemic1, OR epidemic2,...]
d_epiORsd = {} # d_epiORsd[T] = sd of ORs across simulations that resulted in epidemics
d_episize = defaultdict(list) # epidemic sizes of epidemics only
d_ressize = defaultdict(list) # epidemics sizes of all simulations


### parameters ###
numsims = 2000 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
# Tlist = np.linspace(0, .2, num=21, endpoint=True) # probability of transmission
Tlist = np.linspace(0.04, 0.07, num = 16, endpoint=True) # zoom in on realistic values


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
print "number of child and adult nodes:", c_size, a_size


###############################################
### T simulations ###
for T in Tlist:
	print "T value for current simulations:", T
	# d_binlist[simnumber] = [list of 0s and 1s in node numbered index - 1 if node was infected in simnumber simulation]
	d_binlist = defaultdict(list) 
	start = clock()
	for num in np.arange(numsims):
		child_rec, adult_rec, total_rec, bin_list = perc.perc_age_gen(G, d_node_age, T)
		d_simresults[(T, num)] = (child_rec, adult_rec, total_rec)
		d_binlist[num] = bin_list
	print "simtime, simnum:", clock()-start, "\t", num

	# print binary file of infecteds for each set of T simulations
	filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/binlist_T_%ssims_T%s_vax0.txt' %(numsims, str(T))
	pp.print_dictlist_to_file(d_binlist, filename)

##############################################
### calculate OR for all results ###
# # calculate incidence of children and adults
# for key in d_simresults:
# 	d_resincid[(key[0], 'C')].append(d_simresults[key][0]/c_size)
# 	d_resincid[(key[0], 'A')].append(d_simresults[key][1]/a_size)
#
# ### clean up all size results (all results) ###
# for T in Tlist:
# 	d_ressize[T] = [d_simresults[key][2] for key in d_simresults if T == key[0]]
#
# ### calculate OR of children to adults incidence for all epidemics ###
# # calculate OR
# for T in Tlist:
# # 	print T, d_epiincid[(T, 'C')], d_epiincid[(T,'A')]
# 	numerator = [(item/(1-item)) for item in d_resincid[(T, 'C')]]
# 	denominator = [(item/(1-item)) for item in d_resincid[(T, 'A')]]
# 	d_resOR[T] = [n/d for n, d in zip(numerator, denominator) if n != 0 and d != 0]
# 	d_resORsd[T] = np.std(d_resOR[T]) # for errorbars in plot
# 	d_ressize[T] = [d_simresults[key][2] for key in d_simresults if T == key[0]]
# 	# need to remove item from d_ressize list, indexes where numerator or denominator equal 0 because OR cannot be divided by 0. 8/14/13 ECL

# NOTE: how to treat data when item == 1??

##############################################
### calculate OR for epidemics only ###
# # separate epidemics from all results
d_simepi = perc.epidemicsonly(d_simresults, size_epi) 

# calculate incidence of children and adults (epidemics only)
for key in d_simepi.keys():
	d_epiincid[(key[0], 'C')].append(d_simepi[key][0]/c_size)
	d_epiincid[(key[0], 'A')].append(d_simepi[key][1]/a_size)

### calculate OR of children to adults incidence for all epidemics ###
# grab unique list of T values that produced at least one epidemic
T_epi = list(set([key[0] for key in d_simepi]))
# calculate OR
for T in T_epi:
# 	print T, d_epiincid[(T, 'C')], d_epiincid[(T,'A')]
	numerator = [(item/(1-item)) for item in d_epiincid[(T, 'C')]]
	denominator = [(item/(1-item)) for item in d_epiincid[(T, 'A')]]
	d_epiOR[T] = [n/d for n,d in zip(numerator, denominator)]
	d_epiORsd[T] = np.std(d_epiOR[T]) # for errorbars in plot

### clean up episize (epidemics only) ###
T_epi = list(set([key[0] for key in d_simepi]))
for T in T_epi:
	d_episize[T] = [d_simepi[key][2] for key in d_simepi if T == key[0]]
d_episize_sd = dict((k, np.std(d_episize[k])) for k in d_episize)



##############################################
### plot OR by T ###
plt.errorbar(T_epi, [np.mean(d_epiOR[T]) for T in T_epi], yerr=[d_epiORsd[T] for T in T_epi], marker='o', color='black', linestyle='None')
plt.xlabel('T')
plt.ylabel('OR, child:adult')
plt.xlim([-.05, 0.25])
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
plt.savefig(figname)
# plt.show()
plt.close()

# ECL 8/14/13 dictionaries need to be fixed
# ### plot OR by epidemic size (all results) ###
# ressize, resOR = [],[]
# for T in sorted(Tlist):
# 	for item in d_ressize[T]:
# 		ressize.append(item)
# 	for item in d_epiOR[T]:
# 		resOR.append(item)
# plt.scatter(ressize, resOR, marker = 'o', facecolors='none', edgecolors='black', s = 45, label='all sims')
# plt.xlabel('epidemic size')
# plt.ylabel('OR, child:adult')
# plt.legend(loc = 'upper left')
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/resOR_size_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
# # plt.savefig(figname)
# plt.show()
# plt.close()

### plot OR by episize (epidemics only) ###
episize, epiOR = [],[]
for T in sorted(T_epi):
	for item in d_episize[T]:
		episize.append(item)
	for item in d_epiOR[T]:
		epiOR.append(item)
plt.scatter(episize, epiOR, marker = 'o', facecolors='none', edgecolors='blue', s = 45, label='epidemics only')
plt.xlabel('epidemic size')
plt.ylabel('OR, child:adult')
plt.legend(loc = 'upper left')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_size_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
plt.savefig(figname)
# plt.show()
plt.close()

##############################################
### DIAGNOSTICS: epidemic size w/ error bars ###
# plot episize by T
plt.errorbar(T_epi, [np.mean(d_episize[T]) for T in T_epi], yerr=[d_episize_sd[T] for T in T_epi], marker='o', color='black', linestyle='None')
plt.xlabel('T')
plt.ylabel('epidemic size')
plt.xlim([-.05, 0.25])
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/episize_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
plt.savefig(figname)
# plt.show()
plt.close()

### DIAGNOSTICS: number of epidemics ### 
d_numepi ={}
T_epi = list(set([key[0] for key in d_simepi]))
for T in T_epi:
	d_numepi[T] = len([d_simepi[key][2] for key in d_simepi if T == key[0]])

# plot
plt.plot(sorted(T_epi), [d_numepi[T] for T in sorted(T_epi)], marker='o', color='black')
plt.xlabel('T')
plt.ylabel('number of epidemics')
plt.xlim([-.05, 0.25])
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/numepi_T_%ssims_T%s-%s_vax0.png' %(numsims, str(min(Tlist)), str(max(Tlist)))
plt.savefig(figname)
# plt.show()
plt.close()


##############################################
### write dictionaries to files ###
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_T_%ssims_T%s-%s_vax0.txt'%(numsims, str(min(Tlist)), str(max(Tlist)))
pp.print_OR_to_file(d_epiOR, filename)

