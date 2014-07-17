#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 7/16/14
###Function:
##### Visualize results of time-based epidemic simulations for immunity_single sims when aligned by epidemic time, which is defined as aligning tsteps at which simulation attained 5% of cumulative infections during the epidemic.


###Import data: 

###Command Line: python age_time_immunity_single_epitime_viz.py
##############################################

### notes ###
# Ages:
# 1 = Infant, 2 = Toddler, 3 = Child, 4 = Adult, 5 = Senior, 6 = Elder (in nursing home)

# Places (edge attribute):
# F = household/family, S = school, H = hospital, M = shopping mall, W = workplace, D = daycare, E = elsehwere, P = preschool, O = nursing homes, N = neighbor

# T_critical = 0.0565868

### packages/modules ###
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import zipfile
from time import clock
import bisect
import pylab as P

## local modules ##
import percolations as perc
import simulation_parameters as par
import pretty_print as pp

### data processing parameters ###
align_prop = par.dp_alignprop

### plotting parameters ###
numsims = par.pp_numsims
size_epi = par.pp_size_epi
inf_period = par.pp_inf_period
g = par.pp_gamma
T = par.pp_T
b = par.pp_b
# mean retro OR params
beg_perc, end_perc = par.pp_beg_retro_perc, par.pp_end_retro_perc
# specific to immunity params
imm_val = par.pp_immune_val
prop_ls = par.pp_prop_list
zstring = par.pp_pstr_range
zstring2 = par.pp_mstr_fixed

print "Params:", numsims, size_epi, inf_period, g, T, b, imm_val, prop_ls

### data structures ###
d_node_age = {} # d_node_age[nodenumber] = ageclass

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/immunity_time_%ssims_beta%.3f_%s_%s.zip' %(numsims, b, zstring, zstring2)

#############################################
# age data processing
graph_ages = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/urban_network_added_with_info_May24_2014/urban_ages_N10k_Sept2012.txt') # node number and age class

for line in graph_ages:
    new_line = line.strip().split(' ')
    node, age = new_line
    d_node_age[node] = age # node-ageclass dictionary

# define network size
N = len(d_node_age)
print "network size:", N

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
# dict_epiincid[(code, simnumber, 'T', 'C' or 'A')] = [T, C or A incid at tstep 0, T, C or A incid at tstep 1...], where incidence is simply number of new cases (raw)
# dict_epiAR[(code, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per population size
# dict_epiOR[(code, simnumber)] = [OR at tstep0, OR at tstep1...]
# dict_epiOR_filt[(code, simnum)] = [OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers]
# dict_epiresults[(code, simnumber)] = (episize, c_episize, a_episize)
# d_totepiOR[code] = [attack rate OR at sim1, OR at sim 2...]
# d_retroOR[code] = [mean retro OR at sim1, mean retro OR at sim2...]
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt, d_totepiOR, d_retroOR = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)

for prop in prop_ls:
	zstring3 = 'prop%s' %(prop) # string for filename disambiguation
	processing = clock()
	Itstep_file = 'Results/Itstep_immunity_time_%ssims_beta%.3f_%s_%s.txt' %(numsims, b, zstring3, zstring2)
	Rtstep_file = 'Results/Rtstep_immunity_time_%ssims_beta%.3f_%s_%s.txt' %(numsims, b, zstring3, zstring2)
	# recreate epidata from zip archive
	d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata2(Itstep_file, Rtstep_file, zipname, prop, size_epi, ch, ad, to, sr, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
	# calculate OR over entire simulation
	d_totepiOR[prop] = perc.OR_sim(numsims, d_epiresults, prop, chsz, adsz)
	# calculate mean retro OR
	d_retroOR[prop] = perc.mean_retro_OR(d_epiincid, d_epiOR, prop, beg_perc, end_perc)

	print prop, "processed", clock() - processing

	# number of simulations that reached epidemic size
	num_epi = sum([1 for key in d_epiresults if d_epiresults[key][0] > size_epi])
	print prop, "number of epidemics", num_epi
# grab unique list of proportion values that produced at least one epidemic
prop_epi = list(set([key[0] for key in d_epiincid]))

print d_retroOR

#############################################
## draw plots

### plot total simulation AR with SD bars for children, adults, toddlers and the elderly vs pre-existing immunity
c_mns, c_sds, a_mns, a_sds = [],[],[],[]
d_mns, d_sds, s_mns, s_sds = [],[],[],[]

for v in sorted(prop_epi):
	# attack rate by age group
	C_episz_allsims = [sum(d_epiincid[key])/chsz for key in d_epiincid if key[0] == v and key[2] == 'C']
	A_episz_allsims = [sum(d_epiincid[key])/adsz for key in d_epiincid if key[0] == v and key[2] == 'A']
	D_episz_allsims = [sum(d_epiincid[key])/tosz for key in d_epiincid if key[0] == v and key[2] == 'D']
	S_episz_allsims = [sum(d_epiincid[key])/srsz for key in d_epiincid if key[0] == v and key[2] == 'S']
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
CH = plt.errorbar(sorted(prop_epi), c_mns, yerr = c_sds, marker = 'o', color = 'red', linestyle = 'None')
AD = plt.errorbar(sorted(prop_epi), a_mns, yerr = a_sds, marker = 'o', color = 'blue', linestyle = 'None')
TO = plt.errorbar(sorted(prop_epi), d_mns, yerr = d_sds, marker = 'o', color = 'orange', linestyle = 'None')
SR = plt.errorbar(sorted(prop_epi), s_mns, yerr = s_sds, marker = 'o', color = 'green', linestyle = 'None')
plt.xlabel('proportion of adults with any immunity (epidemics only)')
plt.ylabel('Attack Rate')
lines = [CH, AD, TO, SR]
plt.legend(lines, ['children (5-18)', 'adults (19-64)', 'toddlers (0-2)', 'seniors (65+)'], loc = 'upper left')
plt.xlim([0, 1])
plt.ylim([0, 1])
figname = 'Figures/HR-AR_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring, zstring2)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

#############################################
## plot total simulation OR with std bars vs pre-existing immunity

plt.errorbar(sorted(prop_epi), [np.mean(d_totepiOR[val]) for val in sorted(prop_epi)], yerr = [np.std(d_totepiOR[val]) for val in sorted(prop_epi)], marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('proportion of adults with any immunity (epidemics only)')
plt.ylabel(par.pf_OR_lab)
plt.xlim([0, 1])
figname = 'Figures/totepiOR_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring, zstring2)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

#############################################
## plot mean retro OR with std bars vs pre-existing immunity

plt.errorbar(sorted(prop_epi), [np.mean(d_retroOR[val]) for val in sorted(prop_epi)], yerr = [np.std(d_retroOR[val]) for val in sorted(prop_epi)], marker = 'o', color = 'black', linestyle = 'None')
plt.xlabel('proportion of adults with any immunity (epidemics only)')
plt.ylabel(par.pf_mnretro_lab)
plt.xlim([0, 1])
figname = 'Figures/mnretro-epiOR_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring, zstring2)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
for prop in prop_epi:
	zstring3 = 'prop%s' %(prop) # string for filename disambiguation

	##############################################
	### plot filtered and aligned OR by time###
	# alignment at tstep where sim reaches 5% of total episize
	# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
	# each sim is one line, each susc is a diff color on one plot

	ORonly = clock()

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[s] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, prop, align_prop)

	# realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
	avg_align_tstep = 0

	# plot aligned data
	# zip beta, episim number, and tstep for 5% cum-inf for sims where (s, episim number) is the key for d_epiOR_filt
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[prop]):
		plt.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
		plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('epidemic time step, 5-95% cum infections')
	plt.ylabel(par.pf_OR_lab)

	figname = 'Figures/epiORalign_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring3, zstring2)
	plt.savefig(figname)
	plt.clf()
	pp.compress_to_ziparchive(zipname, figname)
	print "ORonly plotting time", prop, clock() - ORonly
	# plt.show()

	##############################################
	### plot filtered and aligned OR by time ###
	### secondary axis with child and adult incidence ###
	# alignment at tstep where sim reaches 5% of total episize
	# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
	# each sim is one line, each beta is a diff color on one plot

	ORincid = clock()

	# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
	# d_dummyalign_tstep[suscept_val] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, prop, align_prop)

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
	for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[prop]):

		## OR y-axis
		OR, = yax_OR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')

		## AR y-axis
		child, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'C')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'C')][t5:]], marker = 'None', color = 'red')
		adult, = yax_AR.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiAR[(k0, k1, 'A')][t5:])), [AR * 100 for AR in d_epiAR[(k0, k1, 'A')][t5:]], marker = 'None', color = 'blue')

	# plot settings
	lines = [OR, child, adult]
	yax_OR.legend(lines, par.pf_epiORincid_leg, loc = 'upper right')
	yax_OR.set_ylabel(par.pf_OR_lab)
	yax_OR.set_xlabel('epidemic time step, 5-95% cum infections')
	yax_AR.set_ylabel('Incidence per 100')

	# save plot
	plt.xlabel('time step')
	figname = 'Figures/epiORincid_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring3, zstring2)
	plt.savefig(figname)
	plt.clf()
	pp.compress_to_ziparchive(zipname, figname)
	print "ORincid plotting time", prop, clock() - ORincid
	# plt.show()











