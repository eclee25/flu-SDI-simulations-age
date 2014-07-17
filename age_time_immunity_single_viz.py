#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 7/16/14
###Function:
##### Visualize results of time-based epidemic simulations where pre-existing immunity exists and is heterogeneous within the adult population. Vary single average value of immunity for subpopulation of adults with any pre-existing immunity.

###Import data: 

###Command Line: python age_time_immunity_single_viz.py
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

## local modules ##
import percolations as perc
import simulation_parameters as par
import pretty_print as pp

### plotting parameters ###
numsims = par.pp_numsims
size_epi = par.pp_size_epi
inf_period = par.pp_inf_period
g = par.pp_gamma
T = par.pp_T
b = par.pp_b
# specific to immunity params
imm_val_ls = par.pp_immune_val_list
prop = par.pp_prop
zstring = par.pp_pstr_fixed
zstring2 = par.pp_mstr_range

print "Params:", numsims, size_epi, inf_period, g, T, b, imm_val_ls, prop

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

##############################################
# data processing - convert tstep info into dictionaries

# declare dictionaries
# dict_epiincid[(code, simnumber, 'T', 'C' or 'A')] = [T, C or A incid at tstep 0, T, C or A incid at tstep 1...], where incidence is simply number of new cases (raw)
# dict_epiAR[(code, simnumber, 'T', 'C' or 'A')] = [T, C or A attack rate at tstep 0, T, C or A attack rate at tstep 1...], where attack rate is number of new cases per population size
# dict_epiOR[(code, simnumber)] = [OR at tstep0, OR at tstep1...]
# dict_epiOR_filt[(code, simnum)] = [OR for each time step for epidemics only where OR is nan when we want to exclude the time point due to small infected numbers]
# dict_epiresults[(code, simnumber)] = (episize, c_episize, a_episize)
d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = defaultdict(list), defaultdict(list), {}, defaultdict(list), defaultdict(list)

for imm_val in imm_val_ls:
	zstring3 = 'single%s' %(imm_val) # string for filename disambiguation
	processing = clock()
	Itstep_file = 'Results/Itstep_immunity_time_%ssims_beta%.3f_%s_%s.txt' %(numsims, b, zstring, zstring3)
	Rtstep_file = 'Results/Rtstep_immunity_time_%ssims_beta%.3f_%s_%s.txt' %(numsims, b, zstring, zstring3)
	d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt = perc.recreate_epidata(Itstep_file, Rtstep_file, zipname, imm_val, size_epi, ch, ad, d_epiincid, d_epiOR, d_epiresults, d_epiAR, d_epiOR_filt)
	print "processed", clock() - processing

	# number of simulations that reached epidemic size
	num_epi = sum([1 for key in d_epiresults if d_epiresults[key][0] > size_epi])
	print imm_val, "number of epidemics", num_epi
	# grab unique list of immunity magnitude values that produced at least one epidemic
	imm_epi = list(set([key[0] for key in d_epiincid]))

#############################################
## draw plots 
for imm in imm_epi:
	pl_ls = [key for key in d_epiOR if key[0] == imm]
	zstring3 = 'single%s' %(imm) # string for filename disambiguation
	##############################################
	### plot OR by time ###
	# each epidemic sim is one line
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR[key])), d_epiOR[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('time step')
	plt.ylabel(par.pf_OR_lab)
	figname = 'Figures/epiOR_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring, zstring3)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
	# plt.show()

for imm in imm_epi:
	pl_ls = [key for key in d_epiOR if key[0] == imm]
	zstring3 = 'single%s' %(imm) # string for filename disambiguation
	##############################################
	### plot filtered OR by time ###
	# each sim is one line
	for key in pl_ls:
		plt.plot(xrange(len(d_epiOR_filt[key])), d_epiOR_filt[key], marker = 'None', color = 'grey')
	plt.plot(xrange(250), [1] * 250, marker = 'None', color = 'red', linewidth = 2)
	plt.xlabel('sim time step, 5-95% cum infections')
	plt.ylabel(par.pf_OR_lab)
	figname = 'Figures/epiORfilt_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring, zstring3)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
	# plt.show()

for imm in imm_epi:
	pl_ls = [key for key in d_epiOR if key[0] == imm]
	zstring3 = 'single%s' %(imm) # string for filename disambiguation
	##############################################
	### plot incidence by time ###
	# each sim is one line
	pl_ls = [key for key in d_epiincid if key[0] == imm and key[2] == 'T']
	for key in pl_ls:
		plt.plot(xrange(len(d_epiincid[key])), d_epiincid[key], marker = 'None', color = 'grey')	
	plt.xlabel('time step')
	plt.ylabel('number of new cases')
	figname = 'Figures/epiincid_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring, zstring3)
	plt.savefig(figname)
	plt.close()
	pp.compress_to_ziparchive(zipname, figname)
	# plt.show()


##############################################
### plot hist of episize by immunity value ###
d_episize = defaultdict(list)
for v in imm_epi:
	d_episize[v] = [sum(d_epiincid[key]) for key in d_epiincid if key[0] == v and key[2] == 'T']
plt.errorbar(imm_epi, [np.mean(d_episize[v]) for v in imm_epi], yerr = [np.std(d_episize[v]) for v in imm_epi], marker = 'o', color = 'black', linestyle = 'None')
plt.xlim([-0.1, 1.1])
plt.xlabel('adult immunity')
plt.ylabel('epidemic size')
figname = 'Figures/episize_immunity_time_%ssims_beta%.3f_%s_%s.png' %(numsims, b, zstring, zstring2)
plt.savefig(figname)
plt.close()
pp.compress_to_ziparchive(zipname, figname)
# plt.show()

