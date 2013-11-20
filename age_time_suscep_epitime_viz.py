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
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import zipfile
from time import clock
import bisect

## local modules ##
import percolations as perc
import pretty_print as pp


### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']

### data processing parameters ###
align_prop = 0.05

### simulation parameters ###
numsims = 800  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 1/float(5) # 5 days recovery here
# T = 0.0643 # total epidemic size = 20%
T = 0.075 

# T = beta / (beta + gamma)
# when T = 0.0643 and gamma = 1/5, b = 0.0137
# when T = 0.075 and gamma = 1/5, b = 0.0162
b = (-T * gamma)/(T - 1) 

# define different child susceptibilities
# Cauchemez 2004 cites child susceptibility to be 1.15 times greater than that of adults
s1, s2 = 1, 1.5
susc_list = np.linspace(s1, s2, num=6, endpoint=True)

### data structures ###
# d_node_age[nodenumber] = ageclass
d_node_age = {} 

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/suscep_time_%ssims_beta%.3f_suscep%.1f-%.1f_vax0.zip' %(numsims, b, s1, s2)

#############################################
# age data processing
graph_ages = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

for line in graph_ages:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age # node-ageclass dictionary

# define network size
N = len(d_node_age)

# create binary lists to indicate children and adults
ch = [1 if d_node_age[str(node)] == '3' else 0 for node in xrange(1, int(N) + 1)]
ad = [1 if d_node_age[str(node)] == '4' else 0 for node in xrange(1, int(N) + 1)]

##############################################
# data processing - convert tstep info into dictionaries

# declare dictionaries
# dict_epiincid[(s, simnumber, 'T', 'C' or 'A')] = [T, C or A incid at tstep 0, T, C or A incid at tstep 1...], where incidence is simply number of new cases (raw)
# dict_epiAR[(s, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per population size
# dict_epiOR[(s, simnumber)] = [OR at tstep0, OR at tstep1...]
# dict_epiOR_filt[(s, simnum)] = [OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers]
# dict_epiresults[(s, simnumber)] = (episize, c_episize, a_episize)
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list)

for s in susc_list:
	processing = clock()
	# reference filenames in zipfolder
	Itstep_file = 'Results/Itstep_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	Rtstep_file = 'Results/Rtstep_susc_time_%ssims_beta%.3f_susc%.1f_vax0.txt' %(numsims, b, s)
	# recreate epidata from zip archive
	d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata(Itstep_file, Rtstep_file, zipname, s, size_epi, ch, ad, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
	print s, "processed", clock() - processing

# grab unique list of child susceptibility values that produced at least one epidemic
susc_epi = list(set([key[0] for key in d_epiincid]))

##############################################
### plot filtered and aligned OR by time for each suscep value ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each susc is a diff color on one plot

for s in susc_epi:
	ORonly = clock()

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[s] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, s, align_prop)

	# TEST (11/19/13): realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
	avg_align_tstep = 0

	# plot aligned data
	# zip beta, episim number, and tstep for 5% cum-inf for sims where (s, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[s]):
	
		plt.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('epidemic time step, adult suscep: ' + str(s) + ', 5-95% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 20])
	plt.xlim([-1, 150])
	figname = 'Figures/epiORalign_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot filtered and aligned OR by time for each suscep value ###
### secondary axis with child and adult incidence ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each beta is a diff color on one plot
 
for s in susc_epi:
	ORincid = clock()

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[suscept_val] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, s, align_prop)

	# TEST (11/19/13): realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
	avg_align_tstep = 0

	# PROCESS YAX_AR: 
	# call upon d_epiAR dictionary
	# dict_epiAR[(s, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per 100 individuals

	# plot data
	# create two y-axes
	fig, yax_OR = plt.subplots()
	yax_AR = yax_OR.twinx()
	
	# zip s, episim number, and tstep for 5% cum-inf for sims where (beta, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[s]):
	
		## OR y-axis
		OR, = yax_OR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
		
		## AR y-axis
		child, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'C')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'C')][t5:]], marker = 'None', color = 'red')
		adult, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'A')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'A')][t5:]], marker = 'None', color = 'blue')

	# plot settings
	lines = [OR, child, adult]
	yax_OR.legend(lines, ['Odds Ratio', 'Child Incidence', 'Adult Incidence'], loc = 'upper right')
	yax_OR.set_ylabel('OR, child:adult')
	yax_OR.set_ylim([0, 20])
	yax_OR.set_xlim([-1, 150])
	yax_OR.set_xlabel('epidemic time step, adult suscep: ' + str(s) + ', 5-95% cum infections')
	yax_AR.set_ylabel('Incidence per 100')
	yax_AR.set_ylim([0, 4])
	
	# save plot
	figname = 'Figures/epiORincid_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()


























