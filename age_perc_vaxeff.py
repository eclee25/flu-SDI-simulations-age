#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/13/13
###Function:
##### 1) odds ratio for two set vax coverage levels (0.245 and 0.455) and susceptibility values that range from 0 to 1 by 0.05, where 0 is fully protected and 1 is fully susceptible. Sweet spot is around 0.1 to 0.5?

###Import data: urban_edges_Sarah.csv, urban_ages_Sarah.csv

###Command Line: python age_perc_vaxeff.py
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

### set simulation parameters ###
numsims = 4000 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
T = 0.065 # ~20% AR in naive population
# Tcov = 0.077 # ~40% AR in naive population
susceplist = np.linspace(0, 1, num=21, endpoint=True) # susceptibility multipliers
vaxcov_scenarios = ['Low', 'High']
# just for reference
vaxcov_vals = [0.245, 0.455]

### import data ###
f = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_edges_Sarah.csv') # Vancouver network
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

### build network ###
G=nx.Graph()
for edge in f:
	G.add_edge(*edge.strip().split(','))
net_size = float(G.order())
print "network size:", net_size

### construct node-age class dictionary ###
for line in file2:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age # node-ageclass dictionary

# number of nodes of each age class
c_size, a_size = perc.child_adult_size(d_node_age)
print "number of child and adult nodes:", c_size, a_size


##############################################
### vax eff simulations ###
for scenario in vaxcov_scenarios:
	vc_scen = scenario
	start = clock()
	
	# simulations
	for s in susceplist:
		print "simulations on susceptibility level: ", s
		# d_binlist[simnumber] = [list of 0s and 1s in node numbered index - 1 if node was infected in simnumber simulation]
		d_binlist = defaultdict(list)
		
		for num in np.arange(numsims):
			child_rec, adult_rec, tot_rec, tot_vax, bin_list = perc.perc_age_vaxeff(G, d_node_age, T, vc_scen, s) 
			d_simresults[(vc_scen, s, num)] = (child_rec, adult_rec, tot_rec, tot_vax)
			d_binlist[num] = bin_list
#		print "simtime, simnum:", clock()-start, "\t", num
	# print binary file of infecteds for each set of cov simulations
	filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_vaxeff_%ssims_T%s_cov%s_suscep%s.txt'%(numsims, str(T), scenario, str(s))
	pp.print_dictlist_to_file(d_binlist, filename)

	
##############################################
### calculate incidence for children and adults for each simulation that turned into an epidemic ###
# separate epidemics from all results
# key = (vc_scen, s, num), val = (child_rec, adult_rec, tot_rec, tot_vax)
d_simepi = perc.epidemicsonly(d_simresults, size_epi) 

# calculate incidence of children and adults
for key in d_simepi:
	# pull dictionary entries
	vc_scen, s = key[0], key[1]
	child_rec, adult_rec = d_simepi[key][0], d_simepi[key][1]
	# calculate incidence
	d_epiincid[(vc_scen, s, 'C')].append(child_rec/c_size)
	d_epiincid[(vc_scen, s, 'A')].append(adult_rec/a_size)

##############################################
### calculate OR of children to adults incidence for all epidemics ###
# grab unique list of vaxcov scenarios and susceptibility combinations that produced at least one epidemic
vc_s_epi = list(set([(key[0], key[1]) for key in d_simepi]))

# calculate OR
for key in vc_s_epi:
	vc_scen, s = key[0], key[1]
	numerator = [(item/(1-item)) for item in d_epiincid[(vc_scen, s, 'C')]]
	denominator = [(item/(1-item)) for item in d_epiincid[(vc_scen, s, 'A')]]
	# ORs of epis where key = vaxcov_scen, s
	d_epiOR[key] = [n/d for n,d in zip(numerator, denominator)]
	print "number of epidemics at", key, ": ", len(d_epiOR[key])
	# dict for errorbars in plot where key = vaxcov_scen, s
	d_epiORsd[key] = np.std(d_epiOR[key]) 

# NOTE: how to treat data when item == 1??

##############################################
### OUT: plot OR by vaxeff ###
# separate lines for different vaxcov scenarios
colorlist = ['black', 'blue']
for scenario in vaxcov_scenarios:
	s_epi = list(set([key[1] for key in d_simepi if scenario == key[0]]))
		
	# plot
	plt.errorbar(s_epi, [np.mean(d_epiOR[(scenario, s)]) for s in s_epi], yerr=[d_epiORsd[(scenario, s)] for s in s_epi], marker='o', color=colorlist.pop(0), linestyle='None', label = scenario + ' coverage')
plt.xlabel('susceptibility (vax eff = 1 - susceptibility)')
plt.ylabel('OR, child:adult')
plt.xlim([-.05, 1.05])
plt.legend(loc = 'upper left')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_vaxeff_%ssims_T%s_cov%s-%s_suscep%s-%s.png'%(numsims, str(T), '0.245', '0.455', str(min(susceplist)), str(max(susceplist)))
plt.savefig(figname)
plt.close()


##############################################
### DIAGNOSTICS: epidemic size w/ error bars ###
colorlist = ['black', 'blue']
for scenario in vaxcov_scenarios:
	d_episize = defaultdict(list)
	s_epi = list(set([key[1] for key in d_simepi if scenario == key[0]]))
	for s in s_epi:
		d_episize[s] = [d_simepi[key][2] for key in d_simepi if scenario == key[0] and s == key[1]]
	d_episize_sd = dict((k, np.std(d_episize[k])) for k in d_episize)
	print scenario, 'episize', d_episize.items()
	# plot
	plt.errorbar(s_epi, [np.mean(d_episize[s]) for s in s_epi], yerr=[d_episize_sd[s] for s in s_epi], marker='o', color=colorlist.pop(0), linestyle='None', label = scenario + ' coverage')
plt.xlabel('susceptibility (vax eff = 1 - susceptibility)')
plt.ylabel('epidemic size')
plt.xlim([-.05, 1.05])
plt.legend(loc = 'upper left')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/diagnostics/episize_vaxeff_%ssims_T%s_cov%s-%s_suscep%s-%s.png'%(numsims, str(T), '0.245', '0.455', str(min(susceplist)), str(max(susceplist)))
plt.savefig(figname)
plt.close()
# plt.show()

### DIAGNOSTICS: number of epidemics ### 
colorlist = ['black', 'blue']
for scenario in vaxcov_scenarios:
	d_numepi = {}
	s_epi = list(set([key[1] for key in d_simepi if scenario == key[0]]))
	for s in s_epi:
		d_numepi[s] = len([d_simepi[key][2] for key in d_simepi if scenario == key[0] and s == key[1]])
		# plot
	plt.plot(sorted(s_epi), [d_numepi[s] for s in sorted(s_epi)], marker='o', color=colorlist.pop(0), label = scenario + ' coverage')
plt.xlabel('susceptibility (vax eff = 1 - susceptibility)')
plt.ylabel('number of epidemics')
plt.xlim([-.05, 1.05])
plt.legend(loc = 'upper left')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/diagnostics/numepi_vaxeff_%ssims_T%s_cov%s-%s_suscep%s-%s.png'%(numsims, str(T), '0.245', '0.455', str(min(susceplist)), str(max(susceplist)))
plt.savefig(figname)
plt.close()
# plt.show()

##############################################
### write dictionaries to files ###
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_vaxeff_%ssims_T%s_cov%s-%s_suscep%s-%s.txt'%(numsims, str(T), str(min(vaxcov_vals)), str(max(vaxcov_vals)), str(min(susceplist)), str(max(susceplist)))
pp.print_OR_to_file(d_epiOR, filename)




