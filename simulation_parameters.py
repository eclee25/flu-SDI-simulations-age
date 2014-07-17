#!/usr/bin/python

##############################################
###Python template
### Author: Elizabeth Lee
### Date: 6/11/14
### List simulation parameters and plotting parameters. Separate parameters are preferable because simulations and plotting may be occurring simultaneously for different sets of experiments.

##############################################
### Packages ###
import numpy as np

### Simulation Parameters ###

## base parameters ##
sp_numsims = 800 # number of simulations
sp_size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)

## disease parameters ##
sp_inf_period = 5. # on avg, assume 5 day infectious period
sp_gamma = 1/sp_inf_period # gamma = probability of recovery at each time step
sp_T = 0.0643 # total epidemic size = 20%
# sp_T = 0.075 # total epidemic size = 30%

# dependent on sp_T and sp_gamma, defined above
sp_b = (-sp_T * sp_gamma)/(sp_T - 1) # T = beta / (beta + gamma)
# when T = 0.0643 and gamma = 1/5, b = 0.0137
# when T = 0.075 and gamma = 1/5, b = 0.0162

### immunity sims ###

##  single sims 
sp_prop = 0.1 # proportion of adults randomly assigned any pre-existing immunity - subpopulation proportion is fixed
sp_immune_val_list = np.linspace(0, 1, 11) # magnitude of pre-existing immunity when setting 'prop' adult proportion of population to an average value - subpop immunity magnitude varies in "_single" sims
sp_pstr_fixed = 'prop%s' %(sp_prop)
sp_mstr_range = 'single%s-%s' %(min(sp_immune_val_list), max(sp_immune_val_list))

## prop sims
sp_prop_list = np.linspace(0, 1, 11) # proportion of adults randomly assigned any pre-existing immunity - subpopulation proportion varies in "_prop" sims
sp_immune_val = 0.2 # magnitude of pre-existing immunity when setting 'prop' adult proportion of population to an average value - subpop immunity magnitude is fixed
sp_pstr_range = 'prop%s-%s' %(min(sp_prop_list), max(sp_prop_list))
sp_mstr_fixed = 'single%s' %(sp_immune_val)

## other sims (not yet implemented)
sp_max_immune_val = 1 # maximum pre-existing immunity magnitude when draw from a uniform distribution - max immunity for uniform draw is fixed
sp_max_immune_val = np.linspace(0, 1, 11) # maximum pre-existing immunity magnitude when draw from a uniform distribution - max immunity for uniform draw varies in "_unif" sims

### infectious period (recovery) sims ###
sp_r1 = 3 # minimum infectious period length
sp_r2 = 15 # maximum infectious period length in days
sp_rnum = 9 # number of denominations in between

### sick behavior sims ###
sp_delay = 1 # days passed of infectious period before contacts are cut 
sp_cut = 0.5 # proportion of random school/work contacts cut for children and adults, respectively, to represent the change in behavior
sp_pop = "A" # string of populations that cut contacts ("C" or "A" for children and adults at school and work, respectively)


##############################################
### Data Processing Parameters ###
dp_alignprop = 0.05 # align plots along "epidemic time," which means the time is aligned at the tstep at which the simulation attained 'dp_alignprop' proportion of cumulative infections over the course of the epidemic


##############################################
### Plotting Parameters ### 

## base parameters ##
pp_numsims = 800 # number of simulations
pp_size_epi = 515 # threshold value that designates an epidemic in the network (5% of network)
	
## disease parameters ##
pp_inf_period = 5. # on avg, assume 5 day infectious period
pp_gamma = 1/sp_inf_period # gamma = probability of recovery at each time step
pp_T = 0.0643 # total epidemic size = 20%
# sp_T = 0.075 # total epidemic size = 30%
pp_b = (-sp_T * sp_gamma)/(sp_T - 1) # T = beta / (beta + gamma)
# when T = 0.0643 and gamma = 1/5, b = 0.0137
# when T = 0.075 and gamma = 1/5, b = 0.0162

## mean retro OR calculation parameters ##
# averages were calculated in SDI_Data/explore/scripts/create_fluseverity_figs/data_extraction/cum_incid_classif_periods.py
pp_beg_retro_perc = 0.39
pp_end_retro_perc = 0.45


### immunity plots ###
##  single sims 
pp_prop = 0.1 # proportion of adults randomly assigned any pre-existing immunity - subpopulation proportion is fixed
pp_immune_val_list = np.linspace(0, 1, 11) # magnitude of pre-existing immunity when setting 'prop' adult proportion of population to an average value - subpop immunity magnitude varies in "_single" sims
pp_pstr_fixed = 'prop%s' %(pp_prop)
pp_mstr_range = 'single%s-%s' %(min(pp_immune_val_list), max(pp_immune_val_list))

## prop sims
pp_prop_list = np.linspace(0, 1, 11) # proportion of adults randomly assigned any pre-existing immunity - subpopulation proportion varies in "_prop" sims
pp_immune_val = 0.2 # magnitude of pre-existing immunity when setting 'prop' adult proportion of population to an average value - subpop immunity magnitude is fixed
pp_pstr_range = 'prop%s-%s' %(min(pp_prop_list), max(pp_prop_list))
pp_mstr_fixed = 'single%s' %(pp_immune_val)

## other sims (not yet implemented)
pp_max_immune_val = 1 # maximum pre-existing immunity magnitude when draw from a uniform distribution - max immunity for uniform draw is fixed
pp_max_immune_val = np.linspace(0, 1, 11) # maximum pre-existing immunity magnitude when draw from a uniform distribution - max immunity for uniform draw varies in "_unif" sims


### infectious period (recovery) plots ###
pp_r1 = 3 # minimum infectious period length
pp_r2 = 15 # maximum infectious period length in days
pp_rnum = 9 # number of denominations in between

### sick behavior plots ###
pp_delay = 1 # days passed of infectious period before contacts are cut 
pp_cut = 0.5 # proportion of random school/work contacts cut for children and adults, respectively, to represent the change in behavior
pp_pop = "A" # string of populations that cut contacts ("C" or "A" for children and adults)

##############################################
### Plot Formatting ###

pf_colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']
pf_HRAR_cols = ['red', 'blue', 'orange', 'green']
pf_HRAR_leg = ['children (5-18)', 'adults (19-64)', 'toddlers (0-2)', 'seniors (65+)']
pf_epiORincid_leg = ['Odds Ratio', 'Child Incidence', 'Adult Incidence']
pf_OR_lab = 'OR, child:adult'
pf_mnretro_lab = 'Mean Retrospective OR'