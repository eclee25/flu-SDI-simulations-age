#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/7/13
###Purpose: visualize results of time-based epidemic simulations when aligned by epidemic time, which is defined as aligning tsteps at which simulation attained 5% of cumulative infections during the epidemic
#### pairs with age_time_suscep.py

###Import data: 

###Command Line: python age_time_suscep_epitime_viz.py
##############################################

####### notes #######
### codebook of age class codes
# '1' - Toddlers: 0-2
# '2' - Preschool: 3-4
# '3' - Children: 5-18
# '4' - Adults: 19-64
# '5' - Seniors: 65+ (community)
# '6' - Elders: 65+ (nursing home)
# There are only 94 "elders" in the Vancouver network, and they all reside in one nursing home, so they can be combined with the seniors for analysis purposes (all_elderly).

### packages/modules ###
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pretty_print as pp
from collections import defaultdict
import zipfile
import percolations as perc
import bisect


### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']

### data processing parameters ###
align_prop = 0.05

### pickled data parameters ###
numsims = 800  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T = 0.0643 (total epidemic size = 20%)
T = 0.0643
# T = beta / (beta + gamma)
b = (-T * gamma)/(T - 1) # b = 0.0137
# define different adult susceptibilities
s1, s2 = 0, 1
susc_list = np.linspace(s1, s2, num=11, endpoint=True)

### import pickled data ###
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname2 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiincid_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epipreval_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/suscepi_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_filt_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_tot_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.txt' %(numsims, b, s1, s2)
d_epiOR = pickle.load(open(pname1, "rb"))
d_epiincid = pickle.load(open(pname2, "rb"))
d_epipreval = pickle.load(open(pname3, "rb"))
susc_epi = pickle.load(open(pname4, "rb"))
d_epiOR_filt = pickle.load(open(pname5, "rb"))
d_epiOR_tot = pickle.load(open(pname6, "rb"))

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/suscep_time_%ssims_beta%.3f_suscep%.1f-%.1f_vax0.zip' %(numsims, b, s1, s2)

### process node ages ###
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class
d_node_age={} # d_node_age[nodenumber] = ageclass

for line in file2:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age # node-ageclass dictionary
# create child and adult indicator vectors where index = nodenumber - 1
c_nodes = [1 if d_node_age[str(node)] == '3' else 0 for node in xrange(1, 10305)]
a_nodes = [1 if d_node_age[str(node)] == '4' else 0 for node in xrange(1, 10305)]


##############################################
### plot filtered and aligned OR by time for each suscep value ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each susc is a diff color on one plot

for s in susc_epi:

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[beta] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, s, align_prop)

	# plot aligned data
	# zip beta, episim number, and tstep for 5% cum-inf for sims where (beta, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[s]):
	
		plt.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('epidemic time step, adult suscep: ' + str(s) + ', 10-90% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 60])
	plt.xlim([-1, 200])
	figname = 'Figures/epiORalign_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
# 	plt.savefig(figname)
# 	plt.close()
# 	pp.compress_to_ziparchive(zipname, figname)
	plt.show()

##############################################
### plot filtered and aligned OR by time for each suscep value ###
### secondary axis with child and adult incidence ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each beta is a diff color on one plot
 
for s in susc_epi:

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[suscept_val] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, s, align_prop)
		
	# PROCESS YAX_AR: 
	# recreate child and adult incid from results txt files
	extractfile = 'Results/Itstep_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	# d_incid[(simnumber, ageclass)] = [num_inf_ageclass at t0,
	d_incid = perc.recreate_incid(extractfile, zipname, size_epi, c_nodes, a_nodes)

	# plot data
	# create two y-axes
	fig, yax_OR = plt.subplots()
	yax_AR = yax_OR.twinx()
	
	# zip beta, episim number, and tstep for 5% cum-inf for sims where (beta, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[s]):
	
		## OR y-axis
		OR, = yax_OR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
		
		## AR y-axis
		child, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_incid[(k1, 'C')][t5:])), d_incid[(k1, 'C')][t5:], marker = 'None', color = 'red')
		adult, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_incid[(k1, 'A')][t5:])), d_incid[(k1, 'A')][t5:], marker = 'None', color = 'blue')

	# plot settings
	lines = [OR, child, adult]
	yax_OR.legend(lines, ['Odds Ratio', 'Child Incidence', 'Adult Incidence'], loc = 'upper right')
	yax_OR.set_ylabel('OR, child:adult')
	yax_OR.set_ylim([0, 60])
	yax_OR.set_xlim([-1, 200])
	yax_OR.set_xlabel('epidemic time step, adult suscep: ' + str(s) + ', 10-90% cum infections')
	yax_AR.set_ylabel('Incidence per 100')
	yax_AR.set_ylim([0, 2])
	
	# save plot
	figname = 'Figures/epiORincid_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()


























