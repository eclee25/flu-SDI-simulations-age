#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/8/13
###Function:
##### 1) odds ratio by vax coverage

###Import data: urban_edges_Sarah.csv, urban_ages_Sarah.csv

###Command Line: python age_perc_vaxcov.py
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
d_epiincid = defaultdict(list) # d_epiincid[(T, 'C' or 'A')] = [child/adult incidence epidemic1, child/adult incidence epidemic2,...]
d_epiOR = defaultdict(list) # d_epiOR[T] = [OR epidemic1, OR epidemic2,...]
d_epiORsd = {} # d_epiORsd[T] = sd of ORs across simulations that resulted in epidemics

### parameters ###
numsims = 500 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
vaxcovlist = np.linspace(0, .2, num=21, endpoint=True) # vax coverage
T = 0.065 # ~20% AR in naive population
# T = 0.077 # ~40% AR in naive population
cov_fixed = 0.245 # reference value for OR vs vaxeff plot

### import data ###
f = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_edges_Sarah.csv') # Vancouver network
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

### construct node-age class dictionary ###
G=nx.Graph()
for edge in f:
	G.add_edge(*edge.strip().split(','))
net_size = float(G.order())
print "network size:", net_size

for line in file2:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age # node-ageclass dictionary

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)
print "number of child and adult nodes:", c_size, a_size


##############################################
### vax cov simulations ###
for cov in vaxcovlist:
	print "vax coverage for current sim:", cov
	# d_binlist[simnumber] = [list of 0s and 1s in node numbered index - 1 if node was infected in simnumber simulation]
	d_binlist = defaultdict(list)
	
	# construct graph w/ rnd vax
	VG=G.copy()
	perc.rnd_vax(VG, cov) # random vax at coverage lvl cov

	# simulations
	start = clock()
	for num in np.arange(numsims):
		child_rec, adult_rec, total_rec, bin_list = perc.perc_age_gen(G, d_node_age, T)
		d_simresults[(cov, num)] = (child_rec, adult_rec, total_rec)
		d_binlist[num] = bin_list
	print "simtime, simnum:", clock()-start, "\t", num

	# print binary file of infecteds for each set of cov simulations
	filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/binlist_cov_%ssims_T%s_eff%s.txt'%(numsims, str(T), str(round(cov/cov_fixed, 3)))
	pp.print_dictlist_to_file(d_binlist, filename)

##############################################
### calculate incidence for children and adults for each simulation that turned into an epidemic ###
# # separate epidemics from all results
d_simepi = perc.epidemicsonly(d_simresults, size_epi) 

# calculate incidence of children and adults
for key in d_simepi:
	cov = key[0]
	child_rec, adult_rec = d_simepi[key][0], d_simepi[key][1]
	d_epiincid[(cov, 'C')].append(child_rec/c_size)
	d_epiincid[(cov, 'A')].append(adult_rec/a_size)

##############################################
### calculate OR of children to adults incidence for all epidemics ###
# grab unique list of coverage values that produced at least one epidemic
cov_epi = sorted(list(set([key[0] for key in d_simepi])))
# convert coverage * efficacy to efficacy only with constant coverage
eff_epi = [cov/cov_fixed for cov in sorted(cov_epi)]
# calculate OR
for cov in sorted(cov_epi):
# 	print cov, d_epiincid[(cov, 'C')], d_epiincid[(cov, 'A')]
	numerator = [(item/(1-item)) for item in d_epiincid[(cov, 'C')]]
	denominator = [(item/(1-item)) for item in d_epiincid[(cov, 'A')]]
	d_epiOR[cov] = [n/d for n,d in zip(numerator, denominator)]
	d_epiORsd[cov] = np.std(d_epiOR[cov]) # for errorbars in plot

# NOTE: how to treat data when item == 1??

##############################################
### plot OR by vaxeff ###
plt.errorbar(eff_epi, [np.mean(d_epiOR[cov]) for cov in sorted(cov_epi)], yerr = [d_epiORsd[cov] for cov in sorted(cov_epi)], marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('random vax efficacy (cov = 0.245)')
plt.ylabel('OR, child:adult')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_cov_%ssims_T%s_eff%s-%s.png'%(numsims, str(T), str(round(min(vaxcovlist)/cov_fixed, 3)), str(round(max(vaxcovlist)/cov_fixed, 3)))
plt.savefig(figname)
# plt.show()
plt.close()

##############################################
### DIAGNOSTICS: epidemic size w/ error bars ###
d_episize = defaultdict(list)
cov_epi = sorted(list(set([key[0] for key in d_simepi])))
for cov in sorted(cov_epi):
	d_episize[cov] = [d_simepi[key][2] for key in d_simepi if cov == key[0]]
d_episize_sd = dict((k, np.std(d_episize[k])) for k in d_episize)

### plot episize by cov ###
plt.errorbar(eff_epi, [np.mean(d_episize[cov]) for cov in sorted(cov_epi)], yerr=[d_episize_sd[cov] for cov in sorted(cov_epi)], marker='o', color='black', linestyle='None')
plt.xlabel('random vax efficacy (cov = 0.245)')
plt.ylabel('epidemic size')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/episize_cov_%ssims_T%s_eff%s-%s.png'%(numsims, str(T), str(round(min(vaxcovlist)/cov_fixed, 3)), str(round(max(vaxcovlist)/cov_fixed, 3)))
plt.savefig(figname)
# plt.show()
plt.close()

### DIAGNOSTICS: number of epidemics ### 
d_numepi ={}
cov_epi = sorted(list(set([key[0] for key in d_simepi])))
for cov in sorted(cov_epi):
	d_numepi[cov] = len([d_simepi[key][2] for key in d_simepi if cov == key[0]])

# plot
plt.plot(eff_epi, [d_numepi[cov] for cov in sorted(cov_epi)], marker='o', color='black')
plt.xlabel('random vax efficacy (cov = 0.245)')
plt.ylabel('number of epidemics')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/numepi_cov_%ssims_T%s_eff%s-%s.png'%(numsims, str(T), str(round(min(vaxcovlist)/cov_fixed, 3)), str(round(max(vaxcovlist)/cov_fixed, 3)))
plt.savefig(figname)
# plt.show()
plt.close()

##############################################
### write dictionaries to files ###
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_cov_%ssims_T%s_eff%s-%s.txt'%(numsims, str(T), str(round(min(vaxcovlist)/cov_fixed, 3)), str(round(max(vaxcovlist)/cov_fixed, 3)))
pp.print_OR_to_file(d_epiOR, filename)



