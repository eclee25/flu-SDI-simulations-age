#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/7/13
###Purpose: visualize results of time-based epidemic simulations with varying susceptibility values
#### pairs with age_time_suscep.py

###Import data: 

###Command Line: python age_time_suscep_viz.py
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
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink', 'brown', 'indigo']

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
### plot OR by time for each suscep value ###
# each epidemic sim is one line
for s in susc_epi:
	pl_ls = [key for key in d_epiOR if key[0] == s]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, adult suscep: ' + str(s))
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 30])
	plt.xlim([-1, 200])
	figname = 'Figures/epiOR_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot filtered OR by time for each suscep value ###
# each sim is one line
for s in susc_epi:
	pl_ls = [key for key in d_epiOR_filt if key[0] == s]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('sim time step, adult suscep: ' + str(s) + ', 5-95% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 30])
	plt.xlim([-1, 200])
	figname = 'Figures/epiORfilt_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot filtered OR by time for all s values ###
# each sim is one line, each suscep is a diff color on one plot
for s in susc_epi:
	pl_ls = [key for key in d_epiOR_filt if key[0] == s]
	colvec = colorvec.pop()
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = colvec)
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, all s values, 5-95% cum infections')
	plt.ylabel('filtered OR, child:adult')
	plt.ylim([0, 30])
	plt.xlim([-1, 200])
figname = 'Figures/epiORfilt_susc_time_%ssims_beta%.3f_allsuscep_vax0.png' %(numsims, b)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot incidence by time for each s value ###
# each sim is one line
for s in susc_epi:
	pl_ls = [key for key in d_epiincid if key[0] == s and key[2] == 'T']
	for key in pl_ls:
		plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')
	plt.xlabel('time step, adult suscep: ' + str(s))
	plt.ylabel('number of new cases')
	plt.xlim([-1, 200])
	figname = 'Figures/epiincid_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot hist of episize by child suscep value ###
d_episize = defaultdict(list)
for s in susc_epi:
	d_episize[s] = [sum(d_epiincid[key]) for key in d_epiincid if key[0] == s and key[2] == 'T']
plt.errorbar(susc_epi, [np.mean(d_episize[s]) for s in susc_epi], yerr = [np.std(d_episize[s]) for s in susc_epi], marker = 'o', color = 'black', linestyle = 'None')
plt.xlim([0.9, 1.6])
plt.xlabel('adult susceptibility')
plt.ylabel('epidemic size')
figname = 'Figures/episize_susc_time_%ssims_beta%.3f_vax0.png' %(numsims, b)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()



# ##############################################
# ### this dict is no longer created in data processing
# ### plot total OR by suscep ###
# # one chart for all adult suscep values
# plt.errorbar(susc_epi, [np.mean(d_epiOR_tot[s]) for s in susc_epi], yerr = [np.std(d_epiOR_tot[s]) for s in susc_epi], marker = 'o', color = 'black', linestyle = 'None')
# plt.xlabel('adult susceptibility')
# plt.ylabel('OR, child:adult')
# plt.ylim([-1, 30])
# figname = 'Figures/epiORtot_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.png' %(numsims, b, s1, s2)
# plt.savefig(figname)
# plt.close()
# pp.compress_to_ziparchive(zipname, figname)
# # plt.show()














