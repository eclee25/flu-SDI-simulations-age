#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/9/13
###Function: 
## plot epidemic sizes and error bars to see amount of variation in sizes and determine if more than 1000 simulations are needed
###Import data: simresults_1000sims_T0.0-0.2_vax0.txt, simresults_1000sims_T0.065_vaxcov0.0-0.2.txt

###Command Line: python 
##############################################


### notes ###


### packages/modules ###
import csv
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

## local modules ##

### data structures ###
d_T_episize = defaultdict(list)
d_cov_episize = defaultdict(list)

### parameters ###

### functions ###

### import data ###
f_T = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/T_1000_T0.0-0.2_vax0/simresults_1000sims_T0.0-0.2_vax0.txt', 'r')
lines_T = f_T.readlines()
f_cov = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Results/cov_4000_T0.065_vaxcov0.0-0.2/simresults_4000sims_T0.065_vaxcov0.0-0.2.txt', 'r')
lines_cov = f_cov.readlines()
### program ###

# T plot, remove print dictionary to file formatting #
# for line in lines_T:
# 	# clean up formatting in file
# 	line_c1 = line.replace('),(', ', ')
# 	line_c2 = line_c1.replace('(', '')
# 	line_c3 = line_c2.replace(')', '')
# 	fields = line_c3.split(', ')
# 	
# 	# input only epidemics into dictionary
# 	if float(fields[4]) > 515:
# 		d_T_episize[fields[0]].append(float(fields[4]))
#
# # sd for epidemics only
# d_T_episize_sd = dict((k, np.std(d_T_episize[k])) for k in d_T_episize)
#
# # list of T values that produced at least one epidemic
# T_epi = list(set([key for key in d_T_episize]))
#
# # plot episize by T
# plt.errorbar(T_epi, [np.mean(d_T_episize[T]) for T in T_epi], yerr = [d_T_episize_sd[T] for T in T_epi], color = 'black', linestyle = 'None', marker = 'o')
# plt.xlabel('T')
# plt.ylabel('epidemic size')
# plt.xlim([0, 0.25])
# plt.show()

########################################################

# vaxcov plot, remove print dictionary to file formatting #
for line in lines_cov:
	# clean up formatting in file
	line_c1 = line.replace('),(', ', ')
	line_c2 = line_c1.replace('(', '')
	line_c3 = line_c2.replace(')', '')
	fields = line_c3.split(', ')
	
	# input only epidemics into dictionary
	if float(fields[4]) > 515:
		d_cov_episize[fields[0]].append(float(fields[4]))




# sd for epidemics only
d_cov_episize_sd = dict((k, np.std(d_cov_episize[k])) for k in d_cov_episize)

# list of cov values that produced at least one epidemic
cov_epi = [key for key in d_cov_episize]

# print d_cov_episize[cov_epi[4]]
# print cov_epi[4], np.mean(d_cov_episize[cov_epi[4]])
# print np.std(d_cov_episize[cov_epi[4]])

# plot episize by vax coverage
plt.errorbar(cov_epi, [np.mean(d_cov_episize[cov]) for cov in cov_epi], yerr = [d_cov_episize_sd[cov] for cov in cov_epi], color = 'black', linestyle = 'None', marker = 'o')
plt.xlabel('proportion of vax coverage')
plt.ylabel('epidemic size')
plt.xlim([0, 0.25])
plt.show()

