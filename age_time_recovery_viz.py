#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 1/1/13
###Function:
##### Visualize results of time-based epidemic simulations where recovery rates vary for children 

###Import data: 

###Command Line: python age_time_recovery_viz.py
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
# T_critical = 0.0565868

### packages/modules ###
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import zipfile
from time import clock

## local modules ##
import percolations as perc
import pretty_print as pp

### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']

### simulation parameters ###
numsims = 800  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
inf_period = 5. # 5 days recovery for all non-children
gamma = 1/inf_period
# T = 0.0643 # total epidemic size = 20%
T = 0.075 # total epidemic size = 30%

# T = beta / (beta + gamma)
# when T = 0.0643 and gamma = 1/5, b = 0.0137
# when T = 0.075 and gamma = 1/5, b = 0.0162
b = (-T * gamma)/(T - 1) 

# define different child recovery rates from 3 to 15 days
r1, r2 = 3, 15
rec_list = np.linspace(r1, r2, num=9, endpoint=True) 

### data structures ###
# d_node_age[nodenumber] = ageclass
d_node_age = {} 

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/recov_time_%ssims_beta%.3f_rec%.1f-%.1f_vax0.zip' %(numsims, b, r1, r2)

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

for r in rec_list:
	processing = clock()
	Itstep_file = 'Results/Itstep_recov_time_%ssims_beta%.3f_rec%.1f_vax0.txt' %(numsims, b, r)
	Rtstep_file = 'Results/Rtstep_recov_time_%ssims_beta%.3f_rec%.1f_vax0.txt' %(numsims, b, r)
	# recreate epidata from zip archive
	d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata(Itstep_file, Rtstep_file, zipname, r, size_epi, ch, ad, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
	print r, "processed", clock() - processing

# grab unique list of child recovery rates that produced at least one epidemic
recov_epi = list(set([key[0] for key in d_epiincid]))

##############################################
### plot OR by time for each child recovery value ###
# each epidemic sim is one line
for r in recov_epi:
	pl_ls = [key for key in d_epiOR if key[0] == r]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, child recovery: ' + str(r))
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 30])
	plt.xlim([-1, 200])
	figname = 'Figures/epiOR_recov_time_%ssims_beta%.3f_recov%.1f_vax0.png' %(numsims, b, r)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot filtered OR by time for each child recovery value ###
# each sim is one line
for r in recov_epi:
	pl_ls = [key for key in d_epiOR_filt if key[0] == r]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('sim time step, child recovery: ' + str(r) + ', 5-95% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 30])
	plt.xlim([-1, 200])
	figname = 'Figures/epiORfilt_recov_time_%ssims_beta%.3f_recov%.1f_vax0.png' %(numsims, b, r)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot filtered OR by time for all recov values ###
# each sim is one line, each suscep is a diff color on one plot
for r in recov_epi:
	pl_ls = [key for key in d_epiOR_filt if key[0] == r]
	colvec = colorvec.pop()
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = colvec)
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, all r values, 5-95% cum infections')
	plt.ylabel('filtered OR, child:adult')
	plt.ylim([0, 30])
	plt.xlim([-1, 200])
figname = 'Figures/epiORfilt_recov_time_%ssims_beta%.3f_allrecov_vax0.png' %(numsims, b)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot incidence by time for each r value ###
# each sim is one line
for r in recov_epi:
	pl_ls = [key for key in d_epiincid if key[0] == r and key[2] == 'T']
	for key in pl_ls:
		plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')	
	plt.xlabel('time step, child recov: ' + str(r))
	plt.ylabel('number of new cases')
	plt.xlim([-1, 200])
	figname = 'Figures/epiincid_recov_time_%ssims_beta%.3f_recov%.1f_vax0.png' %(numsims, b, r)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot hist of episize by child recov value ###
d_episize = defaultdict(list)
for r in recov_epi:
	d_episize[r] = [sum(d_epiincid[key]) for key in d_epiincid if key[0] == r and key[2] == 'T']
plt.errorbar(recov_epi, [np.mean(d_episize[r]) for r in recov_epi], yerr = [np.std(d_episize[r]) for r in recov_epi], marker = 'o', color = 'black', linestyle = 'None')
plt.xlim([2, 16])
plt.xlabel('child recovery')
plt.ylabel('epidemic size')
figname = 'Figures/episize_recov_time_%ssims_beta%.3f_vax0.png' %(numsims, b)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

