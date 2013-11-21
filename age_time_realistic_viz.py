#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 11/21/13
###Function: visualize results of time-based epidemic simulations with realistic child and non-child parameters for susceptibility, infectivity and recovery (sigma, beta, and gamma).

###Import data: urban_ages_Sarah.csv, Itstep_realistic_time..., Rtstep_realistic_time...

###Command Line: python 
##############################################


### notes ###


### packages/modules ###
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import zipfile
from time import clock

## local modules ##
import percolations as perc
import pretty_print as pp

### plotting settings ###

### data structures ###
# d_node_age[nodenumber] = ageclass
d_node_age = {}
# d_age_params[str(age class code)] = (sigma, beta, gamma)
d_age_params = {}
# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else 0]
d_save_I_tstep = defaultdict(list) 
d_save_R_tstep = defaultdict(list)

### simulation parameters ###
numsims = 1000  # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)

## realistic disease parameters (Cauchemez 2004)
# gamma = probability of recovery at each time step
# 3.6 day infectious period for children
# 3.9 day infectious period for non-children
gamma_c, gamma_nc = 1/3.6, 1/3.9  
# sigma = susceptibility of individual
sigma_c, sigma_nc = 1.15, 1.0
# infectivity = probability of infecting a contact over entire infectious period
infectivity_c, infectivity_nc = 0.48, 0.26
infectivity_cnc = infectivity_c/infectivity_nc

## T defined by epidemic size on network
# T = transmissibility of individual
T_nc = 0.0643 # total epidemic size = 20%
# T_nc = 0.075 # total epidemic size = 30%
T_c = T_nc * infectivity_cnc
# T = beta / (beta + gamma)
b_c = (-T_c * gamma_c)/(T_c - 1) 
b_nc = (-T_nc * gamma_nc)/(T_nc - 1)

# assign parameters to dict for ease of importing into a function
c_param_tuple = (sigma_c, b_c, gamma_c)
nc_param_tuple = (sigma_nc, b_nc, gamma_nc)
for ageclass in range(1, 7):
	d_age_params[str(ageclass)] = nc_param_tuple
d_age_params['3'] = c_param_tuple

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.zip' %(numsims, b_nc)

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
# dict_epiincid[(simnumber, 'T', 'C' or 'A')] = [T, C or A incid at tstep 0, T, C or A incid at tstep 1...], where incidence is simply number of new cases (raw)
# dict_epiAR[(simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per population size
# dict_epiOR[(simnumber)] = [OR at tstep0, OR at tstep1...]
# dict_epiOR_filt[(simnum)] = [OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers]
# dict_epiresults[(simnumber)] = (episize, c_episize, a_episize)
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list)

processing = clock()
Itstep_file = 'Results/Itstep_realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.txt' %(numsims, b_nc)
Rtstep_file = 'Results/Rtstep_realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.txt' %(numsims, b_nc)

# recreate epidata from zip archive
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata(Itstep_file, Rtstep_file, zipname, b_nc, size_epi, ch, ad, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
print b_nc, "processed", clock() - processing

##############################################
### plot OR by time ###
# each epidemic sim is one line
for key in d_epiOR:
	plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
plt.plot(xrange(200), [1] * 200, marker = 'None', color = 'red', linewidth = 2)
plt.xlabel('time step, non-child beta: ' + str(round(b_nc, 4)))
plt.ylabel('OR, child:adult')
plt.xlim([0, 100])
figname = 'Figures/epiOR_realistic_time_%ssims_ncbeta%.3f_ncbase.png' %(numsims, b_nc)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot filtered OR by time ###
# each sim is one line
for key in d_epiOR_filt:
	plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = 'grey')
plt.plot(xrange(200), [1] * 200, marker = 'None', color = 'red', linewidth = 2)
plt.xlabel('sim time step, non-child beta: ' + str(round(b_nc, 4)) + ', 5-95% cum infections')
plt.ylabel('OR, child:adult')
plt.xlim([0, 100])
figname = 'Figures/epiORfilt_realistic_time_%ssims_ncbeta%.3f_ncbase.png' %(numsims, b_nc)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot incidence by time ###
# each sim is one line
epiincid_T_keys = [key for key in d_epiincid if key[2] == 'T']
for key in epiincid_T_keys:
	plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')
plt.xlabel('time step,non-child beta: ' + str(round(b_nc, 4)))
plt.ylabel('number of new cases')
plt.xlim([0, 100])

figname = 'Figures/epiincid_realistic_time_%ssims_ncbeta%.3f_ncbase.png' %(numsims, b_nc)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot hist of episize ###
episizes = [d_epiresults[key][0] for key in d_epiresults]
plt.hist(episizes, 10)
plt.ylabel('number of epidemic sims')
plt.xlabel('total attack rate, non-child beta: ' + str(round(b_nc, 4)))
figname = 'Figures/episizehist_realistic_time_%ssims_ncbeta%.3f_ncbase.png' %(numsims, b_nc)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()






























