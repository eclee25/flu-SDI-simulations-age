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
gamma = par.pp_gamma
T = par.pp_T
b = par.pp_b
delay = par.pp_delay
cut = par.pp_cut
pop = par.pp_pop

print "Params:", numsims, size_epi, inf_period, gamma, T, b, delay, cut, pop

### data structures ###
d_node_age = {} # d_node_age[nodenumber] = ageclass
code = str(delay)+str(cut)+str(pop)

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.zip' %(numsims, b, delay, cut, pop)

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
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt, d_totepiOR = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list), defaultdict(list)

processing = clock()
Itstep_file = 'Results/Itstep_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.txt' %(numsims, b, delay, cut, pop)
Rtstep_file = 'Results/Rtstep_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.txt' %(numsims, b, delay, cut, pop)
# recreate epidata from zip archive
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata2(Itstep_file, Rtstep_file, zipname, code, size_epi, ch, ad, to, sr, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
# calculate OR over entire simulation
d_totepiOR[code] = perc.OR_sim(numsims, d_epiresults, code, chsz, adsz)
print code, "processed", clock() - processing

# number of simulations that reached epidemic size
print "number of epidemics", len(d_epiincid)

##############################################
### plot total simulation AR with SD bars for children, adults, toddlers and the elderly vs infectious period (1/gamma)

# attack rate by age group
C_episz_allsims = [sum(d_epiincid[key])/chsz for key in d_epiincid if key[0] == code and key[2] == 'C']
A_episz_allsims = [sum(d_epiincid[key])/adsz for key in d_epiincid if key[0] == code and key[2] == 'A']
D_episz_allsims = [sum(d_epiincid[key])/tosz for key in d_epiincid if key[0] == code and key[2] == 'D']
S_episz_allsims = [sum(d_epiincid[key])/srsz for key in d_epiincid if key[0] == code and key[2] == 'S']

P.figure()
n, bins, patches = P.hist([C_episz_allsims, A_episz_allsims, D_episz_allsims, S_episz_allsims], 10, histtype='bar', color=par.pf_HRAR_cols, label=par.pf_HRAR_leg)
P.legend()
P.xlim([0, 0.5])
P.xlabel('Attack Rate')
P.ylabel('Frequency')
figname = 'Figures/HR-AR_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.png' %(numsims, b, delay, cut, pop)
P.savefig(figname)
P.clf()
pp.compress_to_ziparchive(zipname, figname)
# P.show()

#############################################
## plot total simulation OR hist

plt.hist(d_totepiOR[code])
plt.xlabel(par.pf_OR_lab)
plt.ylabel('frequency')
figname = 'Figures/totepiOR_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.png' %(numsims, b, delay, cut, pop)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

##############################################
### plot filtered and aligned OR by time###
# alignment at tstep where sim reaches 5% of total episize
# starting tstep on plot is mode of tsteps where sim reaches 5% of total episize
# each sim is one line, each susc is a diff color on one plot

ORonly = clock()

# PROCESS X-AXIS: identify tstep at which sim reaches 5% of cum infections for the epidemic
# d_dummyalign_tstep[s] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, code, align_prop)

# realign plots for epitime to start at t = 0 by reassigning avg_align_tstep
avg_align_tstep = 0

# plot aligned data
# zip beta, episim number, and tstep for 5% cum-inf for sims where (s, episim number) is the key for d_epiOR_filt
for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[code]):

	plt.plot(xrange(avg_align_tstep, avg_align_tstep+len(d_epiOR_filt[(k0, k1)][t5:])), d_epiOR_filt[(k0, k1)][t5:], marker = 'None', color = 'grey')
plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
plt.xlabel('epidemic time step, 5-95% cum infections')
plt.ylabel(par.pf_OR_lab)

figname = 'Figures/epiORalign_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.png' %(numsims, b, delay, cut, pop)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
print "ORonly plotting time", code, clock() - ORonly
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
d_dummyalign_tstep, avg_align_tstep, dummyk =  perc.define_epi_time(d_epiincid, code, align_prop)

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
for k0, k1, t5 in zip((k[0] for k in dummyk), (k[1] for k in dummyk), d_dummyalign_tstep[code]):

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
figname = 'Figures/epiORincid_sickbehav_time_%ssims_beta%.3f_delay%s_cut%.2f_%spop.png' %(numsims, b, delay, cut, pop)
plt.savefig(figname)
plt.clf()
pp.compress_to_ziparchive(zipname, figname)
print "ORincid plotting time", code, clock() - ORincid
# plt.show()












