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
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pretty_print as pp
from collections import defaultdict
import zipfile
import percolations as perc


### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']


### pickled data parameters ###
numsims = 1000 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T ranges from 0.0 to 0.2, gamma = 1/5 and T = beta / (beta + gamma)
b1, b2 = (-0.0 * gamma)/(0 - 1), (-0.2 * gamma)/(0.2 - 1) # 0, .05
blist = np.linspace(b1, b2, num=11, endpoint=True) # probability of transmission

### import pickled data ###
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname2 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiincid_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epipreval_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/betaepi_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_filt_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_tot_beta_time_%ssims_beta%.3f-%.3f_vax0' %(numsims, b1, b2)
d_epiOR = pickle.load(open(pname1, "rb"))
d_epiincid = pickle.load(open(pname2, "rb"))
d_epipreval = pickle.load(open(pname3, "rb"))
beta_epi = pickle.load(open(pname4, "rb"))
d_epiOR_filt = pickle.load(open(pname5, "rb"))
d_epiOR_tot = pickle.load(open(pname6, "rb"))

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/beta_time_%ssims_beta%.3f-beta%.3f_vax0.zip' %(numsims, b1, b2)

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
	plt.ylim([-3, 15])
	plt.xlim([-1, 125])
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
	plt.xlabel('sim time step, beta: ' + str(beta) + ', 10-90% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 5])
	plt.xlim([-1, 125])
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
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step, all betas, 10-90% cum infections')
	plt.ylabel('filtered OR, child:adult')
	plt.ylim([0, 5])
	plt.xlim([-1, 150])
figname = 'Figures/epiORfilt_beta_time_%ssims_allbetas_vax0.png' %(numsims)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)

# plt.show()


##############################################
### plot total OR by beta ###
# one chart for all betas (comparable to OR vs T)
plt.errorbar(beta_epi, [np.mean(d_epiOR_tot[b]) for b in beta_epi], yerr = [np.std(d_epiOR_tot[b]) for b in beta_epi], marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('beta')
plt.ylabel('OR, child:adult')
plt.xlim([0.01, 0.06])
figname = 'Figures/epiORtot_beta_time_%ssims_beta%.3f-%.3f_vax0.png' %(numsims, b1, b2)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()


##############################################
### plot incidence by time for each beta value ###
# each sim is one line
for beta in beta_epi:
	pl_ls = [key for key in d_epiincid if key[0] == beta]
	for key in pl_ls:
		plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')
	plt.xlabel('time step, beta: ' + str(beta))
	plt.ylabel('number of new cases')
	plt.xlim([-1, 125])
	figname = 'Figures/epiincid_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()

##############################################
### plot prevalence by time for each beta value ###
# each sim is one line
for beta in beta_epi:
	pl_ls = [key for key in d_epipreval if key[0] == beta]
	for key in pl_ls:
		plt.plot(xrange(len(d_epipreval[key])), d_epipreval[key], marker = 'None', color = 'grey')
	plt.xlabel('time step, beta: ' + str(beta))
	plt.ylabel('total number of cases')
	plt.xlim([-1, 125])
	figname = 'Figures/epipreval_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()
















