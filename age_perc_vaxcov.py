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
numsims = 4000 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
vaxcovlist = np.linspace(0, .2, num=21, endpoint=True) # vax coverage
Tcov = 0.065 # fixed T value for vax coverage simulations


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
	
	# construct graph w/ rnd vax
	VG=G.copy()
	perc.rnd_vax(VG, cov) # random vax at coverage lvl cov

	# simulations
	for num in np.arange(numsims):
		start = clock()
		d_simresults[(cov, num)] = perc.perc_age_gen(VG, d_node_age, Tcov) # value is tuple: (child_rec, adult_rec, total_rec)
# 		print "simtime, simnum:", clock()-start, "\t", num


### calculate incidence for children and adults for each simulation that turned into an epidemic ###
# # separate epidemics from all results
d_simepi = perc.epidemicsonly(d_simresults, size_epi) 

# calculate incidence of children and adults
for key in d_simepi.keys():
	d_epiincid[(key[0], 'C')].append(d_simepi[key][0]/c_size)
	d_epiincid[(key[0], 'A')].append(d_simepi[key][1]/a_size)


### calculate OR of children to adults incidence for all epidemics ###
# grab unique list of coverage values that produced at least one epidemic
cov_epi = list(set([key[0] for key in d_simepi.keys()]))
# calculate OR
for cov in cov_epi:
# 	print cov, d_epiincid[(cov, 'C')], d_epiincid[(cov, 'A')]
	numerator = [(item/(1-item)) for item in d_epiincid[(cov, 'C')]]
	denominator = [(item/(1-item)) for item in d_epiincid[(cov, 'A')]]
	d_epiOR[cov] = [n/d for n,d in zip(numerator, denominator)]
	d_epiORsd[cov] = np.std(d_epiOR[cov]) # for errorbars in plot

# NOTE: how to treat data when item == 1??

### write dictionaries to files ###
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/simresults_%ssims_T%s_vaxcov%s-%s.txt'%(numsims, str(Tcov), str(min(vaxcovlist)), str(max(vaxcovlist)))
pp.print_dictionary_to_file(d_simresults, filename)
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/simepi_%ssims_T%s_vaxcov%s-%s.txt'%(numsims, str(Tcov), str(min(vaxcovlist)), str(max(vaxcovlist)))
pp.print_dictionary_to_file(d_simepi, filename)
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiincid_%ssims_T%s_vaxcov%s-%s.txt'%(numsims, str(Tcov), str(min(vaxcovlist)), str(max(vaxcovlist)))
pp.print_dictionary_to_file(d_epiincid, filename)
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_%ssims_T%s_vaxcov%s-%s.txt'%(numsims, str(Tcov), str(min(vaxcovlist)), str(max(vaxcovlist)))
pp.print_dictionary_to_file(d_epiOR, filename)

### plot OR by T ###
plt.errorbar(cov_epi, [np.mean(d_epiOR[cov]) for cov in cov_epi], yerr=[d_epiORsd[cov] for cov in cov_epi], marker='o', color='black', linestyle='None')
plt.xlabel('random vax coverage proportion')
plt.ylabel('OR, child:adult')
plt.xlim([-.05, 0.25])
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_%ssims_T%s_vaxcov%s-%s.png'%(numsims, str(Tcov), str(min(vaxcovlist)), str(max(vaxcovlist)))
plt.savefig(figname)







