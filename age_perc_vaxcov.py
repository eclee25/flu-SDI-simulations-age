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
import pickle

## local modules ##
import percolations as perc

###############################################
### data structures ###
d_node_age={} # d_node_age[nodenumber] = ageclass
d_simresults = {} # d_simresults[(T, simnumber)] = (number of infected children, number of infected adults)
d_simOR = {} # d_simOR[(cov, simnumber)] = ORval
d_epiOR = defaultdict(list) # d_epiOR[cov] = [OR epidemic1, OR epidemic2,...]

# for diagnostics
d_episize = defaultdict(list) # d_episize[cov] = [episize1, episize2, etc...]
d_numepi ={} # d_numepi[cov] = number of epidemics produced from that coverage level

###############################################
### parameters ###
numsims = 50 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
vaxcovlist = np.linspace(0, .2, num=21, endpoint=True) # vax coverage
# T = 0.065 # ~20% AR in naive population
T = 0.077 # ~40% AR in naive population
cov_fixed = 0.245 # reference value for OR vs vaxeff plot

###############################################
### import data ###
f = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_edges_Sarah.csv') # Vancouver network
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

###############################################
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
### vax efficacy simulations ###
for cov in vaxcovlist:
	print "vax coverage for current sim:", cov
	# d_binlist[simnumber] = [list of 0s and 1s in node numbered index - 1 if node was infected in simnumber simulation]
	d_binlist = defaultdict(list)
	
	# construct graph w/ rnd vax
	VG=G.copy()
	# random vax at coverage lvl cov
	perc.rnd_vax(VG, cov) 

	# simulations
	start = clock()
	for num in xrange(numsims):
		child_rec, adult_rec, total_rec, bin_list, OR_val = perc.perc_age_gen(VG, d_node_age, T)
		d_simresults[(cov, num)] = (child_rec, adult_rec, total_rec)
		d_binlist[num] = bin_list
		d_simOR[(cov, num)] = OR_val
	print "simtime, simnum:", clock()-start, "\t", num

	# print binary file of infecteds for each set of cov simulations
	filename = '/home/elee/Documents/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/binlist_cov_%ssims_T%.3f_eff%.3f.txt' %(numsims, T, (cov/cov_fixed))
# 	pp.print_dictlist_to_file(d_binlist, filename)

##############################################
### subset: epidemics only ###

# subset epidemics from all results
# d_simepi[(cov, num)] = (child_rec, adult_rec, total_rec)
d_simepi = perc.epidemicsonly(d_simresults, size_epi)

# subset ORs for sims that produced epidemics
# d_epiOR[cov] = [OR epi1, OR epi2,...]
for key in sorted(d_simepi):
	d_epiOR[key[0]].append(d_simOR[key])

##############################################
### clean x-axis values ###

# unique coverage values that produced at least one epidemic
cov_epi = sorted(list(set([key[0] for key in d_simepi])))

# convert coverage * efficacy to efficacy only with constant coverage
eff_epi = [cov/cov_fixed for cov in cov_epi]

##############################################
### write dictionaries to files ###
filename = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/epiOR_cov_%ssims_T%.3f_eff%.3f-%.3f.txt' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
# pp.print_OR_to_file(d_epiOR, filename)

##############################################
### pickle
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/covepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/effepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_episize_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_simepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname7 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_numepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pickle.dump(d_epiOR, open(pname1, "wb"))
pickle.dump(cov_epi, open(pname3, "wb"))
pickle.dump(eff_epi, open(pname4, "wb"))
pickle.dump(d_episize, open(pname5, "wb"))
pickle.dump(d_simepi, open(pname6, "wb"))
pickle.dump(d_numepi, open(pname7, "wb"))

########################################################
# instead of using below code, use age_perc_vaxcov_viz.py module

# ##############################################
# ### RESULTS: OR by vax efficacy w/ error bars ###
#
# # plot
# plt.errorbar(eff_epi, [np.mean(d_epiOR[cov]) for cov in cov_epi], yerr = [np.std(d_epiOR[cov]) for cov in cov_epi], marker = 'o', color = 'black', linestyle = 'None')
# plt.xlabel('random vax efficacy (cov = 0.245)')
# plt.ylabel('OR, child:adult')
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_cov_%ssims_T%.3f_eff%.3f-%.3f.png' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
# plt.savefig(figname)
# plt.close()
# # plt.show()
#
# ##############################################
# ### DIAGNOSTICS: epidemic size w/ error by T ###
#
# # grab episize info from epidemic results
# for cov in sorted(cov_epi):
# 	d_episize[cov] = [d_simepi[key][2] for key in d_simepi if cov == key[0]]
#
# ### plot episize by vax efficacy ###
# plt.errorbar(eff_epi, [np.mean(d_episize[cov]) for cov in cov_epi], yerr=[np.std(d_episize[cov]) for cov in cov_epi], marker='o', color='black', linestyle='None')
# plt.xlabel('random vax efficacy (cov = 0.245)')
# plt.ylabel('epidemic size')
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/episize_cov_%ssims_T%.3f_eff%.3f-%.3f.png'%(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
# plt.savefig(figname)
# plt.close()
# # plt.show()
#
# ##############################################
# ### DIAGNOSTICS: number of epidemics by T ### 
#
# # grab number of epidemics from epidemic results
# for cov in sorted(cov_epi):
# 	d_numepi[cov] = len([d_simepi[key][2] for key in d_simepi if cov == key[0]])
#
# # plot number of epidemics by vax efficacy
# plt.plot(eff_epi, [d_numepi[cov] for cov in cov_epi], marker='o', color='black')
# plt.xlabel('random vax efficacy (cov = 0.245)')
# plt.ylabel('number of epidemics')
# figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/numepi_cov_%ssims_T%.3f_eff%.3f-%.3f.png' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
# plt.savefig(figname)
# plt.close()
# # plt.show()





