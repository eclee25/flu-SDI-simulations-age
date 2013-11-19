#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/23/13
###Purpose: visualize results of time-based epidemic simulations
#### pairs with age_perc_T_time.py

###Import data: 

###Command Line: python age_perc_T_time_viz.py
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
import math
from time import clock


### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']

### data parameters ###
numsims = 1000 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T ranges from 0.0 to 0.2, gamma = 1/5 and T = beta / (beta + gamma)
# T1, T2 = 0.0, 0.2
T1, T2 = 0.075, 0.075
# T1, T2 = 0.0643, 0.0643
b1, b2 = (-T1 * gamma)/(T1 - 1), (-T2 * gamma)/(T2 - 1) # 0, .05
blist = np.linspace(b1, b2, num=1, endpoint=True) # probability of transmission

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
	# recreate epidata from zip archive (5 to 45 sec)
	d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata(Itstep_file, Rtstep_file, zipname, beta, size_epi, ch, ad, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
	print beta, "processed", clock() - processing

# grab unique list of betas that produced at least one epidemic
beta_epi = list(set([key[0] for key in d_epiincid]))


##############################################
### plot OR by time for each beta value ###
# each sim is one line
for beta in beta_epi:
	pl_ls = [key for key in d_epiOR if key[0] == beta]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, beta: ' + str(beta))
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 15])
	plt.xlim([0, 200])
	figname = 'Figures/epiOR_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot filtered OR by time for each beta value ###
# each sim is one line
for beta in beta_epi:
	pl_ls = [key for key in d_epiOR_filt if key[0] == beta]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('sim time step, beta: ' + str(beta) + ', 5-95% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 10])
	plt.xlim([0, 200])
	figname = 'Figures/epiORfilt_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot filtered OR by time for all beta values ###
# each sim is one line, each beta is a diff color on one plot
for beta in beta_epi:
	pl_ls = [key for key in d_epiOR_filt if key[0] == beta]
	colvec = colorvec.pop()
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = colvec)
	plt.plot(xrange(200), [1] * len(xrange(200)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, all betas, 5-95% cum infections')
	plt.ylabel('filtered OR, child:adult')
	plt.ylim([0, 8])
	plt.xlim([0, 140])
figname = 'Figures/epiORfilt_beta_time_%ssims_allbetas_vax0.png' %(numsims)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot incidence by time for each beta value ###
# each sim is one line
for beta in beta_epi:
	pl_ls = [key for key in d_epiincid if key[0] == beta and key[2] == 'T']
	for key in pl_ls:
		plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')
	plt.xlabel('time step, beta: ' + str(beta))
	plt.ylabel('number of new cases')
	plt.xlim([-1, 200])
	figname = 'Figures/epiincid_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot histogram of total episizes for each beta value ###
for beta in beta_epi:
	episizes = [sum(d_epiincid[key]) for key in d_epiincid if key[0] == beta and key[2] == 'T']
	plt.hist(episizes, 10)
	plt.xlim([1000,10000])
	plt.ylim([0, 150])
	plt.ylabel('number of epidemic sims')
	plt.xlabel('total attack rate, beta: ' + str(beta))
	figname = 'Figures/episizehist_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()



# ##############################################
# ### this dict is no longer created in data processing
# ### plot total OR by beta ###
# # one chart for all betas (comparable to OR vs T)
# plt.errorbar(beta_epi, [np.mean(d_epiOR_tot[b]) for b in beta_epi], yerr = [np.std(d_epiOR_tot[b]) for b in beta_epi], marker = 'o', color = 'black', linestyle = 'None')
# plt.xlabel('beta')
# plt.ylabel('OR, child:adult')
# plt.xlim([0.01, 0.06])
# figname = 'Figures/epiORtot_beta_time_%ssims_beta%.3f-%.3f_vax0.png' %(numsims, b1, b2)
# plt.savefig(figname)
# plt.close()
# pp.compress_to_ziparchive(zipname, figname)
# # plt.show()












