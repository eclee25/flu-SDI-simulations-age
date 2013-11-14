#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/3/13
###Purpose: visualize results of time-based epidemic simulations when aligned by epidemic time, which is defined as aligning tsteps at which simulation attained 5% of cumulative infections during the epidemic
#### pairs with age_perc_T_time.py

###Import data: 

###Command Line: python age_perc_T_epitime_viz.py
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
import pretty_print as pp
from collections import defaultdict
import zipfile
import percolations as perc
from time import clock


### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']

### data processing parameters ###
align_prop = 0.05

### pickled data parameters ###
numsims = 1000 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T ranges from 0.0 to 0.2, gamma = 1/5 and T = beta / (beta + gamma)
T1, T2 = 0.0, 0.2
# T1, T2 = 0.075, 0.075
# T1, T2 = 0.0643, 0.0643
b1, b2 = (-T1 * gamma)/(T1 - 1), (-T2 * gamma)/(T2 - 1) # 0, .05
blist = np.linspace(b1, b2, num=11, endpoint=True) # probability of transmission

# data structures
# d_node_age[str(node)] = age class
d_node_age = {}

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/beta_time_%ssims_beta%.3f-%.3f_vax0.zip' %(numsims, b1, b2)

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
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list)
# 	# dict_epiincid[(beta, simnumber, 'T', 'C' or 'A')] = [T, C or A incid at tstep 0, T, C or A incid at tstep 1...], where incidence is simply number of new cases (raw)
# 	dict_epiincid = defaultdict(list)
# 	# dict_epiAR[(beta, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per population size
# 	dict_epiAR = defaultdict(list)
# 	# dict_epiOR[(beta, simnumber)] = [OR at tstep0, OR at tstep1...]
# 	dict_epiOR = defaultdict(list)
# 	# dict_epiOR_filt[(beta, simnum)] = [OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers]
# 	dict_epiOR_filt = defaultdict(list)
# 	# dict_epiresults[(beta, simnumber)] = (episize, c_episize, a_episize)
# 	dict_epiresults = {}

for beta in blist:
	processing = clock()
	# reference filenames in zipfolder
	Itstep_file = 'Results/Itstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	Rtstep_file = 'Results/Rtstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	# recreate epidata from zip archive
	d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata(Itstep_file, Rtstep_file, zipname, beta, size_epi, ch, ad, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
	print beta, "processed", clock() - processing

# grab unique list of betas that produced at least one epidemic
beta_epi = list(set([key[0] for key in d_epiincid]))

##############################################
### plot filtered and aligned OR by time for each beta value ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each beta is a diff color on one plot

for beta in beta_epi:
	ORonly = clock()
	
	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[beta] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, beta, align_prop)

	# plot aligned data
	# zip beta, episim number, and tstep for 5% cum-inf for sims where (beta, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[beta]):
	
		plt.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('epidemic time step, beta: ' + str(beta) + ', 5-95% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 8])
	plt.xlim([-1, 100])
	
	print "OR only", beta, clock() - ORonly

	# save plot
	figname = 'Figures/epiORalign_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()
	

##############################################
### plot filtered and aligned OR by time for each beta value ###
### secondary axis with child and adult incidence ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each beta is a diff color on one plot
 
for beta in beta_epi:
	ORincid = clock()

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[beta] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, beta, align_prop)
		
	# PROCESS YAX_AR: 
	# call upon d_epiAR dictionary
	# dict_epiAR[(beta, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per 100 individuals

	# plot data
	# create two y-axes
	fig, yax_OR = plt.subplots()
	yax_AR = yax_OR.twinx()
	
	# zip beta, episim number, and tstep for 10% cum-inf for sims where (beta, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[beta]):
	
		## OR y-axis
		OR, = yax_OR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
		
		## AR y-axis
		child, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'C')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'C')][t5:]], marker = 'None', color = 'red')
		adult, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'A')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'A')][t5:]], marker = 'None', color = 'blue')

	# plot settings
	lines = [OR, child, adult]
	yax_OR.legend(lines, ['Odds Ratio', 'Child Incidence', 'Adult Incidence'], loc = 'upper right')
	yax_OR.set_ylabel('OR, child:adult')
	yax_OR.set_ylim([0, 8])
	yax_OR.set_xlim([-1, 100])
	yax_OR.set_xlabel('epidemic time step, beta: ' + str(beta) + ', 5-95% cum infections')
	yax_AR.set_ylabel('Incidence per 100')
	yax_AR.set_ylim([0, 8])
	
	print "ORincid", beta, clock() - ORonly
	
	# save plot
	figname = 'Figures/epiORincid_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()


























