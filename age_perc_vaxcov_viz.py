#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/23/13
###Purpose: visualize results of time-based epidemic simulations
#### pairs with age_perc_T_time.py

###Import data: pickled datasets

###Command Line: python age_perc_T_time.py
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

### pickled data parameters ###
numsims = 50 # number of simulations
vaxcovlist = np.linspace(0, .2, num=21, endpoint=True) # vax coverage
# T = 0.065 # ~20% AR in naive population
T = 0.077 # ~40% AR in naive population
cov_fixed = 0.245 # reference value for OR vs vaxeff plot

### import pickled data ###
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/covepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/effepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_episize_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_simepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
pname7 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_numepi_cov_%ssims_T%.3f_eff%.3f-%.3f' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
d_epiOR = pickle.load(open(pname1, "rb"))
cov_epi = pickle.load(open(pname3, "rb"))
eff_epi = pickle.load(open(pname4, "rb"))
d_episize = pickle.load(open(pname5, "rb"))
d_simepi = pickle.load(open(pname6, "rb"))
d_numepi = pickle.load(open(pname7, "rb"))

##############################################
### RESULTS: OR by vax efficacy w/ error bars ###

# plot
plt.errorbar(eff_epi, [np.mean(d_epiOR[cov]) for cov in cov_epi], yerr = [np.std(d_epiOR[cov]) for cov in cov_epi], marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('random vax efficacy (cov = 0.245)')
plt.ylabel('OR, child:adult')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_cov_%ssims_T%.3f_eff%.3f-%.3f.png' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
# plt.savefig(figname)
# plt.close()
plt.show()

##############################################
### DIAGNOSTICS: epidemic size w/ error by T ###

# grab episize info from epidemic results
for cov in sorted(cov_epi):
	d_episize[cov] = [d_simepi[key][2] for key in d_simepi if cov == key[0]]

### plot episize by vax efficacy ###
plt.errorbar(eff_epi, [np.mean(d_episize[cov]) for cov in cov_epi], yerr=[np.std(d_episize[cov]) for cov in cov_epi], marker='o', color='black', linestyle='None')
plt.xlabel('random vax efficacy (cov = 0.245)')
plt.ylabel('epidemic size')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/episize_cov_%ssims_T%.3f_eff%.3f-%.3f.png'%(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
# plt.savefig(figname)
# plt.close()
plt.show()

##############################################
### DIAGNOSTICS: number of epidemics by T ### 

# grab number of epidemics from epidemic results
for cov in sorted(cov_epi):
	d_numepi[cov] = len([d_simepi[key][2] for key in d_simepi if cov == key[0]])

# plot number of epidemics by vax efficacy
plt.plot(eff_epi, [d_numepi[cov] for cov in cov_epi], marker='o', color='black')
plt.xlabel('random vax efficacy (cov = 0.245)')
plt.ylabel('number of epidemics')
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/numepi_cov_%ssims_T%.3f_eff%.3f-%.3f.png' %(numsims, T, (min(vaxcovlist)/cov_fixed), (max(vaxcovlist)/cov_fixed))
# plt.savefig(figname)
# plt.close()
plt.show()



