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

### data structures ###
# d_node_age[nodenumber] = ageclass
d_node_age = {}
# d_age_params[str(age class code)] = (sigma, beta, gamma)
d_age_params = {}
# d_save_I_tstep[simnumber] (or d_save_R_tstep) = [time step of infection/recovery where index = node number - 1 else 0]
d_save_I_tstep = defaultdict(list) 
d_save_R_tstep = defaultdict(list)

### data processing parameters ###
align_prop = 0.05

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
# dict_epiincid[(b_nc, simnumber, 'T', 'C' or 'A')] = [T, C or A incid at tstep 0, T, C or A incid at tstep 1...], where incidence is simply number of new cases (raw)
# dict_epiAR[(b_nc, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per population size
# dict_epiOR[(b_nc, simnumber)] = [OR at tstep0, OR at tstep1...]
# dict_epiOR_filt[(b_nc, simnum)] = [OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers]
# dict_epiresults[(b_nc, simnumber)] = (episize, c_episize, a_episize)
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list)

processing = clock()
Itstep_file = 'Results/Itstep_realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.txt' %(numsims, b_nc)
Rtstep_file = 'Results/Rtstep_realistic_time_%ssims_ncbeta%.3f_ncbase_Cauchemez04.txt' %(numsims, b_nc)

# recreate epidata from zip archive
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata(Itstep_file, Rtstep_file, zipname, b_nc, size_epi, ch, ad, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
print b_nc, "processed", clock() - processing

##############################################
### plot filtered and aligned OR by time ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line

ORonly = clock()
# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
# d_dummyalign_tstep[s] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, b_nc, align_prop)

# TEST: realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
avg_align_tstep = 0

# plot aligned data
# zip beta, episim number, and tstep for 5% cum-inf for sims where (s, episim number) is the key for d_epiOR_filt
for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[b_nc]):
	plt.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
plt.xlabel('epidemic time step, non-child beta: ' + str(round(b_nc, 4)) + ', 5-95% cum infections')
plt.ylabel('OR, child:adult')
plt.ylim([0, 15])
plt.xlim([-1, 50])
figname = 'Figures/epiORalign_realistic_time_%ssims_ncbeta%.3f_ncbase.png' %(numsims, b_nc)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
print "ORonly plotting time", clock() - ORonly
# plt.show()

##############################################
### plot filtered and aligned OR by time for each suscep value ###
### secondary axis with child and adult incidence ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each beta is a diff color on one plot
 
ORincid = clock()

# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
# d_dummyalign_tstep[suscept_val] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, b_nc, align_prop)

# TEST: realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
avg_align_tstep = 0

# PROCESS YAX_AR: 
# call upon d_epiAR dictionary
# dict_epiAR[(b_nc, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per 100 individuals

# plot data
# create two y-axes
fig, yax_OR = plt.subplots()
yax_AR = yax_OR.twinx()
	
# zip s, episim number, and tstep for 5% cum-inf for sims where (b_nc, episim number) is the key for d_epiOR_filt
for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[b_nc]):
	
	## OR y-axis
	OR, = yax_OR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
		
	## AR y-axis
	child, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'C')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'C')][t5:]], marker = 'None', color = 'red')
	adult, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'A')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'A')][t5:]], marker = 'None', color = 'blue')

# plot settings
lines = [OR, child, adult]
yax_OR.legend(lines, ['Odds Ratio', 'Child Incidence', 'Adult Incidence'], loc = 'upper right')
yax_OR.set_ylabel('OR, child:adult')
yax_OR.set_ylim([0, 15])
yax_OR.set_xlim([-1, 50])
yax_OR.set_xlabel('epidemic time step, non-child beta: ' + str(round(b_nc, 4)) + ', 5-95% cum infections')
yax_AR.set_ylabel('Incidence per 100')
yax_AR.set_ylim([0, 8])

# save plot
figname = 'Figures/epiORincid_realistic_time_%ssims_ncbeta%.3f_ncbase.png' %(numsims, b_nc)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
print "ORincid plotting time", clock() - ORonly
# plt.show()




























