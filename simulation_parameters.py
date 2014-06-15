#!/usr/bin/python

##############################################
###Python template
### Author: Elizabeth Lee
### Date: 6/11/14
### List simulation parameters and plotting parameters. Separate parameters are preferable because simulations and plotting may be occurring simultaneously for different sets of experiments.

##############################################
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

### sick behavior sims ###
sp_delay = 1 # days passed of infectious period before contacts are cut 
sp_cut = 0.3 # proportion of random school/work contacts cut for children and adults, respectively, to represent the change in behavior
sp_pop = "C" # string of populations that cut contacts ("C" or "A" for children and adults at school and work, respectively)

### infectious period (recovery) sims ###
sp_r1 = 3 # minimum infectious period length
sp_r2 = 15 # maximum infectious period length in days
sp_rnum = 9 # number of denominations in between


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

### sick behavior plots ###
pp_delay = 1 # days passed of infectious period before contacts are cut 
pp_cut = 0.3 # proportion of random school/work contacts cut for children and adults, respectively, to represent the change in behavior
pp_pop = "C" # string of populations that cut contacts ("C" or "A" for children and adults)

### infectious period (recovery) sims ###
pp_r1 = 3 # minimum infectious period length
pp_r2 = 15 # maximum infectious period length in days
pp_rnum = 9 # number of denominations in between

##############################################
### Plot Formatting ###

pf_colorvec = ['black', 'red', 'orange', 'gold', 'green', 'blue', 'cyan', 'darkviolet', 'hotpink']
pf_HRAR_cols = ['red', 'blue', 'orange', 'green']
pf_HRAR_leg = ['children (5-18)', 'adults (19-64)', 'toddlers (0-2)', 'seniors (65+)']
pf_epiORincid_leg = ['Odds Ratio', 'Child Incidence', 'Adult Incidence']
pf_OR_lab = 'OR, child:adult'