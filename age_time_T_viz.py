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


### plotting settings ###
colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']


### data parameters ###
numsims = 100 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
# gamma = probability of recovery at each time step
# on avg, assume 5 days till recovery
gamma = 0.2
# assume T ranges from 0.0 to 0.2, gamma = 1/5 and T = beta / (beta + gamma)
# T1, T2 = 0.0, 0.2
# T1, T2 = 0.075, 0.075
T1, T2 = 0.0643, 0.0643
b1, b2 = (-T1 * gamma)/(T1 - 1), (-T2 * gamma)/(T2 - 1) # 0, .05
blist = np.linspace(b1, b2, num=1, endpoint=True) # probability of transmission

# data structures
# d_Itstep[(beta, simnumber)] = [infection tstep for node 1, infection tstep for node 2, ...]
d_Itstep = defaultdict(list)
# d_Rtstep[(beta, simnumber)] = [recovery tstep for node 1, recovery tstep for node 2, ...]
d_Rtstep = defaultdict(list)
# d_node_age[str(node)] = age class
d_node_age = {}

### ziparchive to read and write results ###
zipname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/beta_time_%ssims_beta%.3f-%.3f_vax0.zip' %(numsims, b1, b2)

graph_ages = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv') # node number and age class

for line in graph_ages:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age # node-ageclass dictionary


##############################################
# data processing - convert tstep info to these dictionaries
for beta in blist:
	# reference filenames in zipfolder
	Itstep_file = 'Results/Itstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	Rtstep_file = 'Results/Rtstep_beta_time_%ssims_beta%.3f_vax0.txt' %(numsims, beta)
	# open ziparchive
	with zipfile.ZipFile(zipname, 'r') as zf:
		simnumber = 0
		# open files in ziparchive
		Itstep = zf.open(Itstep_file, 'r')
		Rtstep = zf.open(Rtstep_file, 'r')
		for line1, line2 in zip(Itstep, Rtstep):
			print line1
			d_Itstep[(beta, simnumber)] = [int(t) if int(t) else 0 for t in line1.split(', ')]
			d_Rtstep[(beta, simnumber)] = [0 if t == 'nan' else int(t) for t in line2.split(', ')]
			simnumber += 1

print d_Itstep.values()

# ch = [1 if d_node_age[str(node)] == '3' else 0 for node in xrange(1, int(net_size) + 1)]
# ad = [1 if d_node_age[str(node)] == '4' else 0 for node in xrange(1, int(net_size) + 1)]
#
# d_incid, d_OR, d_simresults = recreate_simdata(extractfile, zipname, size_epi, 
#
# d_simresults[(beta, num)] = (child_rec, adult_rec, total_rec)
# 		d_simOR[(beta, num)] = OR_list
# 		d_simincid[(beta, num)] = tot_incid_list
# 		d_simpreval[(beta, num)] = tot_preval_list
# 		d_simOR_filt[(beta, num)] = filt_OR_list
# 		d_simOR_tot[(beta, num)] = OR_tot
#
# ## subset: epidemics only ###
# # # subset epidemics from all results
# # key = (beta, simnumber), value = (child_rec, adult_rec, total_rec)
# d_simepi = perc.epidemicsonly(d_simresults, size_epi)
#
# # subset OR/incidence/prevalence values that produced epidemics
# # key = (beta, simnumber), value = [ORs/new cases/total cases/filtered ORs at different timesteps]
# for key in d_simepi:
# 	d_epiOR[key] = d_simOR[key]
# 	d_epiincid[key] = d_simincid[key]
# 	d_epipreval[key] = d_simpreval[key]
# 	d_epiOR_filt[key] = d_simOR_filt[key]
# 	d_epiOR_tot[key[0]].append(d_simOR_tot[key])
#
# # grab unique list of betas that produced at least one epidemic
# beta_epi = list(set([key[0] for key in d_simepi]))

# ##############################################
# ### plot OR by time for each beta value ###
# # each sim is one line
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epiOR if key[0] == beta]
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
# 	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
# 	plt.xlabel('time step, beta: ' + str(beta))
# 	plt.ylabel('OR, child:adult')
# 	plt.ylim([-3, 15])
# 	plt.xlim([-1, 125])
# 	figname = 'Figures/epiOR_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
# 	plt.savefig(figname)
# 	plt.close()
# 	pp.compress_to_ziparchive(zipname, figname)
# # 	plt.show()
#
# ##############################################
# ### plot filtered OR by time for each beta value ###
# # each sim is one line
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epiOR_filt if key[0] == beta]
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = 'grey')
# 	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
# 	plt.xlabel('sim time step, beta: ' + str(beta) + ', 10-90% cum infections')
# 	plt.ylabel('OR, child:adult')
# 	plt.ylim([0, 5])
# 	plt.xlim([-1, 125])
# 	figname = 'Figures/epiORfilt_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
# 	plt.savefig(figname)
# 	plt.close()
# 	pp.compress_to_ziparchive(zipname, figname)
# # 	plt.show()
#
#
# ##############################################
# ### plot filtered OR by time for all beta values ###
# # each sim is one line, each beta is a diff color on one plot
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epiOR_filt if key[0] == beta]
# 	colvec = colorvec.pop()
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = colvec)
# 	plt.plot(xrange(250), [1] * len(xrange(250)), marker = 'None', color = 'red', linewidth = 2)
# 	plt.xlabel('time step, all betas, 10-90% cum infections')
# 	plt.ylabel('filtered OR, child:adult')
# 	plt.ylim([0, 5])
# 	plt.xlim([-1, 150])
# figname = 'Figures/epiORfilt_beta_time_%ssims_allbetas_vax0.png' %(numsims)
# plt.savefig(figname)
# plt.close()
# pp.compress_to_ziparchive(zipname, figname)
#
# # plt.show()
#
#
# ##############################################
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
#
#
# ##############################################
# ### plot incidence by time for each beta value ###
# # each sim is one line
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epiincid if key[0] == beta]
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')
# 	plt.xlabel('time step, beta: ' + str(beta))
# 	plt.ylabel('number of new cases')
# 	plt.xlim([-1, 125])
# 	figname = 'Figures/epiincid_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
# 	plt.savefig(figname)
# 	plt.close()
# 	pp.compress_to_ziparchive(zipname, figname)
# # 	plt.show()
#
# ##############################################
# ### plot prevalence by time for each beta value ###
# # each sim is one line
# for beta in beta_epi:
# 	pl_ls = [key for key in d_epipreval if key[0] == beta]
# 	for key in pl_ls:
# 		plt.plot(xrange(len(d_epipreval[key])), d_epipreval[key], marker = 'None', color = 'grey')
# 	plt.xlabel('time step, beta: ' + str(beta))
# 	plt.ylabel('total number of cases')
# 	plt.xlim([-1, 125])
# 	figname = 'Figures/epipreval_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
# 	plt.savefig(figname)
# 	plt.close()
# 	pp.compress_to_ziparchive(zipname, figname)
# # 	plt.show()

##############################################
### plot histogram of total episizes for each beta value ###
for beta in beta_epi:
	ARs = [sum(d_epiincid[key]) for key in d_epiincid if key[0] == beta]
	plt.hist(ARs, 10)
	plt.ylabel('number of epidemic sims')
	plt.xlabel('total attack rate, beta: ' + str(beta))
	figname = 'Figures/episizehist_beta_time_%ssims_beta%.3f_vax0.png' %(numsims, beta)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
# 	plt.show()















