#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 12/31/13
###Function:
##### Visualize results of time-based epidemic simulations where recovery rates vary for children when aligned by epidemic time, which is defined as aligning tsteps at which simulation attained 5% of cumulative infections during the epidemic


###Import data: 

###Command Line: python age_time_childrecovery_epitime_viz.py
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
inf_period = 5. # 5 days recovery for all non-children
gamma = 1/inf_period
T = 0.0643 # total epidemic size = 20%
# T = 0.075 # total epidemic size = 30%

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
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/childrecov_time_%ssims_beta%.3f_rec%.1f-%.1f_vax0.zip' %(numsims, b, r1, r2)

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
# child and adult population sizes
chsz = float(sum(ch))
adsz = float(sum(ad))

# high risk groups: toddlers (0-2), seniors & elderly (65+)
to = [1 if d_node_age[str(node)] == '1' else 0 for node in xrange(1, int(N) + 1)]
sr = [1 if d_node_age[str(node)] == '5' or d_node_age[str(node)] == '6' else 0 for node in xrange(1, int(N) + 1)]
tosz = float(sum(to))
srsz = float(sum(sr))

print 'children, adults, toddlers, seniors', chsz, adsz, tosz, srsz

##############################################
# data processing - convert tstep info into dictionaries

# storage dictionaries need to be initialized outside the loop
# dict_epiincid[(r, simnumber, 'T', 'C' or 'A')] = [T, C or A incid at tstep 0, T, C or A incid at tstep 1...], where incidence is simply number of new cases (raw)
# dict_epiAR[(r, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per population size
# dict_epiOR[(r, simnumber)] = [OR at tstep0, OR at tstep1...]
# dict_epiOR_filt[(r, simnum)] = [OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers]
# dict_epiresults[(r, simnumber)] = (episize, c_episize, a_episize)
# d_totepiOR[r] = [attack rate OR at sim1, OR at sim 2...]
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt, d_totepiOR = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list), defaultdict(list)

for r in rec_list:
	processing = clock()
	Itstep_file = 'Results/Itstep_childrecov_time_%ssims_beta%.3f_rec%.1f_vax0.txt' %(numsims, b, r)
	Rtstep_file = 'Results/Rtstep_childrecov_time_%ssims_beta%.3f_rec%.1f_vax0.txt' %(numsims, b, r)
	# recreate epidata from zip archive
	d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata2(Itstep_file, Rtstep_file, zipname, r, size_epi, ch, ad, to, sr, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
	# calculate OR over entire simulation
	d_totepiOR[r] = perc.OR_sim(numsims, d_epiresults, r, chsz, adsz)
	print r, "processed", clock() - processing

# grab unique list of child recovery rates that produced at least one epidemic
recov_epi = list(set([key[0] for key in d_epiincid]))
print recov_epi

##############################################
### plot total simulation AR with SD bars for children, adults, toddlers and the elderly vs infectious period (1/gamma)
c_mns, c_sds, a_mns, a_sds = [],[],[],[]
d_mns, d_sds, s_mns, s_sds = [],[],[],[]

for r in sorted(recov_epi):
	# attack rate by age group
	C_episz_allsims = [sum(d_epiincid[key])/chsz for key in d_epiincid if key[0] == r and key[2] == 'C']
	A_episz_allsims = [sum(d_epiincid[key])/adsz for key in d_epiincid if key[0] == r and key[2] == 'A']
	D_episz_allsims = [sum(d_epiincid[key])/tosz for key in d_epiincid if key[0] == r and key[2] == 'D']
	S_episz_allsims = [sum(d_epiincid[key])/srsz for key in d_epiincid if key[0] == r and key[2] == 'S']
	# add mean and SD attack rates to list for each Tmult value
	c_mns.append(np.mean(C_episz_allsims))
	a_mns.append(np.mean(A_episz_allsims))
	d_mns.append(np.mean(D_episz_allsims))
	s_mns.append(np.mean(S_episz_allsims))
	c_sds.append(np.std(C_episz_allsims))
	a_sds.append(np.std(A_episz_allsims))
	d_sds.append(np.std(D_episz_allsims))
	s_sds.append(np.std(S_episz_allsims))
# plot AR by age group with errorbars
CH = plt.errorbar(sorted(recov_epi), c_mns, yerr = c_sds, marker = 'o', color = 'red', linestyle = 'None')
AD = plt.errorbar(sorted(recov_epi), a_mns, yerr = a_sds, marker = 'o', color = 'blue', linestyle = 'None')
TO = plt.errorbar(sorted(recov_epi), d_mns, yerr = d_sds, marker = 'o', color = 'orange', linestyle = 'None')
SR = plt.errorbar(sorted(recov_epi), s_mns, yerr = s_sds, marker = 'o', color = 'green', linestyle = 'None')
plt.xlabel('child infectious period in days (epidemics only)')
plt.ylabel('Attack Rate')
lines = [CH, AD, TO, SR]
plt.legend(lines, ['children (5-18)', 'adults (19-64)', 'toddlers (0-2)', 'seniors (65+)'], loc = 'upper left')
plt.xlim([3, 15])
plt.ylim([0, 1])
figname = 'Figures/HR-AR_childrecov_time_%ssims_beta%.3f_recov%.1f_vax0.png' %(numsims, b, r)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

#############################################
## plot total simulation OR with std bars vs infectious period

# 3/28 troubleshooting
print "sim OR y-axis", [np.mean(d_totepiOR[r]) for r in sorted(recov_epi)]

plt.errorbar(sorted(recov_epi), [np.mean(d_totepiOR[r]) for r in sorted(recov_epi)], yerr = [np.std(d_totepiOR[r]) for r in sorted(recov_epi)], marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('child infectious period in days (epidemics only)')
plt.ylabel('simulation OR, child:adult')
plt.ylim([0, 20])
plt.xlim([3, 15])
figname = 'Figures/totepiOR_childrecov_time_%ssims_beta%.3f_recov_vax0.png' %(numsims, b)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()



#############################################
## plot avg of ORs at each tstep with std bars vs infectious period
mns, sds = [],[]
for r in sorted(recov_epi):
	mns_allsims = [np.mean(np.ma.masked_array(d_epiOR[key], np.isnan(d_epiOR[key]))) for key in d_epiOR if key[0] == r]
	mns.append(np.mean(mns_allsims))
	sds.append(np.mean(mns_allsims))
plt.errorbar(sorted(recov_epi), mns, yerr = sds, marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('child infectious period in days (epidemics only)')
plt.ylabel('simulation OR (avg of avgs), child:adult')
plt.ylim([0, 5])
plt.xlim([3, 15])
figname = 'Figures/totepiOR-avgs_childrecov_time_%ssims_beta%.3frecov%.1f_vax0.png' %(numsims, b, r)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot filtered and aligned OR by time for each child recovery value ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each susc is a diff color on one plot

for r in recov_epi:
	ORonly = clock()

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[s] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, r, align_prop)

	# realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
	avg_align_tstep = 0

	# plot aligned data
	# zip beta, episim number, and tstep for 5% cum-inf for sims where (s, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[r]):
	
		plt.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('epidemic time step, child recov: ' + str(r) + ', 5-95% cum infections')
	plt.ylabel('OR, child:adult')
	plt.ylim([0, 20])
	plt.xlim([-1, 150])
	figname = 'Figures/epiORalign_childrecov_time_%ssims_beta%.3f_recov%.1f_vax0.png' %(numsims, b, r)
	plt.savefig(figname)
	plt.clf()
	pp.compress_to_ziparchive(zipname, figname)
	print "ORonly plotting time", r, clock() - ORonly
# 	plt.show()

##############################################
### plot filtered and aligned OR by time for each suscep value ###
### secondary axis with child and adult incidence ###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each beta is a diff color on one plot
 
for r in recov_epi:
	ORincid = clock()

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[suscept_val] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, r, align_prop)

	# realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
	avg_align_tstep = 0

	# PROCESS YAX_AR: 
	# call upon d_epiAR dictionary
	# dict_epiAR[(r, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per 100 individuals

	# plot data
	# create two y-axes
	fig, yax_OR = plt.subplots()
	yax_AR = yax_OR.twinx()
	
	# zip s, episim number, and tstep for 5% cum-inf for sims where (s, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[r]):
	
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
	yax_OR.set_xlabel('epidemic time step, child recovery: ' + str(r) + ', 5-95% cum infections')
	yax_AR.set_ylabel('Incidence per 100')
	yax_AR.set_ylim([0, 4])
	
	# save plot
	figname = 'Figures/epiORincid_childrecov_time_%ssims_beta%.3f_recov%.1f_vax0.png' %(numsims, b, r)
	plt.savefig(figname)
	plt.clf()
	pp.compress_to_ziparchive(zipname, figname)
	print "ORincid plotting time", r, clock() - ORincid
# 	plt.show()












