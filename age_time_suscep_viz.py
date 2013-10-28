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

### packages/modules ###
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pretty_print as pp
from collections import defaultdict
import zipfile
import percolations as perc


### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink', 'brown', 'indigo']


### pickled data parameters ###
numsims = 800  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 1/float(3)
T = 0.0643 # total epidemic size = 30%
# T = 0.075 # with larger T, adult susceptibility of 0.6 will still be greater than T_critical
# T = beta / (beta + gamma)
b = (-T * gamma)/(T - 1) # when T = 0.0643, b = 0.0137
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

##############################################
### plot OR by time for each suscep value ###
# each sim is one line
for s in susc_epi:
	pl_ls = [key for key in d_epiOR if key[0] == s]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, adult suscep: ' + str(s))
	plt.ylabel('OR, child:adult')
	plt.ylim([-3, 60])
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
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('sim time step, adult suscep: ' + str(s) + ', 10-90% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 60])
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
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, all s values, 10-90% cum infections')
	plt.ylabel('filtered OR, child:adult')
	plt.ylim([0, 60])
	plt.xlim([-1, 200])
figname = 'Figures/epiORfilt_susc_time_%ssims_beta%.3f_allsuscep_vax0.png' %(numsims, b)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)

# plt.show()


##############################################
### plot total OR by suscep ###
# one chart for all adult suscep values
plt.errorbar(susc_epi, [np.mean(d_epiOR_tot[s]) for s in susc_epi], yerr = [np.std(d_epiOR_tot[s]) for s in susc_epi], marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('adult susceptibility')
plt.ylabel('OR, child:adult')
plt.ylim([-1, 30])
figname = 'Figures/epiORtot_susc_time_%ssims_beta%.3f_susc%.1f-%.1f_vax0.png' %(numsims, b, s1, s2)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()


##############################################
### plot incidence by time for each s value ###
# each sim is one line
for s in susc_epi:
	pl_ls = [key for key in d_epiincid if key[0] == s]
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
### plot prevalence by time for each s value ###
# each sim is one line
for s in susc_epi:
	pl_ls = [key for key in d_epipreval if key[0] == s]
	for key in pl_ls:
		plt.plot(xrange(len(d_epipreval[key])), d_epipreval[key], marker = 'None', color = 'grey')
	plt.xlabel('time step, adult suscep: ' + str(s))
	plt.ylabel('total number of cases')
	plt.xlim([-1, 200])
	figname = 'Figures/epipreval_susc_time_%ssims_beta%.3f_susc%.1f_vax0.png' %(numsims, b, s)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot episize by adult susceptibility
d_episize = defaultdict(list)
for s in susc_epi:
	d_episize[s] = [sum(d_epiincid[key]) for key in d_epiincid if key[0] == s]
plt.errorbar(susc_epi, [np.mean(d_episize[s]) for s in susc_epi], yerr = [np.std(d_episize[s]) for s in susc_epi], marker = 'o', color = 'black', linestyle = 'None')
plt.xlim([-0.1, 1.1])
plt.xlabel('adult susceptibility')
plt.ylabel('epidemic size')
figname = 'Figures/episize_susc_time_%ssims_beta%.3f_vax0.png' %(numsims, b)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()
	















