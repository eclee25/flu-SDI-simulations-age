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
import matplotlib.cm as cm

### pickled data parameters ###
numsims = 800 # number of simulations
size_epi = 515 # threshold value that designates an epidemic in the network
Tlist = np.linspace(0.0643, .075, num=2, endpoint=True) # probability of transmission
# Tlist = np.linspace(0.04, 0.07, num = 16, endpoint=True) # zoom in on realistic values

### import pickled data ###
pname1 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_epiOR_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname2 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_episize_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname3 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_ressize_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname4 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_simOR2_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname5 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/Tlist_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname6 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/Tepi_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname7 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_numepi_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))
pname8 = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Pickled/d_simepi_T_%ssims_T%.3f-%.3f_vax0' %(numsims, min(Tlist), max(Tlist))

d_epiOR = pickle.load(open(pname1, "rb"))
d_episize = pickle.load(open(pname2, "rb"))
d_ressize = pickle.load(open(pname3, "rb"))
d_simOR2 = pickle.load(open(pname4, "rb"))
Tlist = pickle.load(open(pname5, "rb"))
T_epi = pickle.load(open(pname6, "rb"))
d_numepi = pickle.load(open(pname7, "rb"))
d_simepi = pickle.load(open(pname8, "rb"))

##############################################
### RESULTS: OR by T w/ error bars ###

# plot
plt.errorbar(T_epi, [np.mean(d_epiOR[T]) for T in T_epi], yerr=[np.std(d_epiOR[T]) for T in T_epi], marker='o', color='black', linestyle='None')
plt.xlabel('T')
plt.ylabel('OR, child:adult')
plt.xlim([0, 0.25])
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_T_%ssims_T%.3f-%.3f_vax0.png' %(numsims, min(Tlist), max(Tlist))
# plt.savefig(figname)
# plt.close()
plt.show()

##############################################
### RESULTS: OR by episize (epidemics only) ###

# colormap formatting
colors = cm.rainbow(np.linspace(0, 1, len(T_epi)))

# separate information from dicts for each T
for T, col in zip(T_epi, colors):
	episize = []
	# d_episize[T] = [episize1, episize2,...]
	for item in d_episize[T]:
		episize.append(item)

	# one plot per T
	# note: d_epiOR[T] = [OR epidemic1, OR epidemic2,...]
	lab_dummy = 'T=' + str(T)
	plt.scatter(episize, d_epiOR[T], marker = 'o', facecolors = 'none', edgecolors = col, s = 45, label = lab_dummy)
plt.xlabel('epidemic size, epidemics only')
plt.ylabel('OR, child:adult')
plt.xlim([0, 10000])
plt.ylim([0, 10])
plt.legend(loc = 'upper left', prop = {'size':8})
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/epiOR_size_T_%ssims_T%.3f_vax0.png' %(numsims, T)
# plt.savefig(figname)
# plt.close()
plt.show()
	
##############################################
### RESULTS: OR by episize (all results) ###

# colormap formatting
colors = cm.rainbow(np.linspace(0, 1, len(Tlist)))

# separate information from dicts for each T
for T, col in zip(Tlist, colors):
	# one plot per T value
	lab_dummy = 'T=' + str(T)
	plt.scatter(d_ressize[T], d_simOR2[T], marker = 'o', facecolors = 'none', edgecolors = col, s = 45, label = lab_dummy)
plt.xlabel('epidemic size, all sims')
plt.ylabel('OR, child:adult')
plt.xlim([0, 10000])
plt.ylim([0, 10])
plt.legend(loc = 'upper center', prop = {'size':6})
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/simOR_size_T_%ssims_T%.3f_vax0.png' %(numsims, T)
# plt.savefig(figname)
# plt.close()
plt.show()

##############################################
### DIAGNOSTICS: epidemic size w/ error by T ###

# plot episize by T
plt.errorbar(T_epi, [np.mean(d_episize[T]) for T in T_epi], yerr=[np.std(d_episize[T]) for T in T_epi], marker='o', color='black', linestyle='None')
plt.xlabel('T')
plt.ylabel('epidemic size')
plt.xlim([0, 0.25])
plt.ylim([0, 10000])
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/episize_T_%ssims_T%.3f-%.3f_vax0.png' %(numsims, min(Tlist), max(Tlist))
# plt.savefig(figname)
# plt.close()
plt.show()

##############################################
### DIAGNOSTICS: number of epidemics by T ### 
for T in T_epi:
	d_numepi[T] = len([d_simepi[key][2] for key in d_simepi if T == key[0]])

# plot
plt.plot(sorted(T_epi), [d_numepi[T] for T in sorted(T_epi)], marker='o', color='black')
plt.xlabel('T')
plt.ylabel('number of epidemics')
plt.xlim([0, 0.25])
figname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Figures/numepi_T_%ssims_T%.3f-%.3f_vax0.png' %(numsims, min(Tlist), max(Tlist))
# plt.savefig(figname)
# plt.close()
plt.show()