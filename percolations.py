#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/7/13
###Purpose: run epidemic simulations on an age-structured network
#### general percolation simulations
#### percolations with random vaccination
#### time-based epidemic simulations
#### calculation of OR from results

###Import data: 

###Command Line: 
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

### structure of dictionaries
# dict_node_age[node] = 'ageclasscode'
# dict_simresults[(T, simulationnumber)] = (recovered_children, recovered_adults, recovered_total)
# dict_simepi[(T, simulationnumber)] = (recovered_children, recovered_adults, recovered_total) # subset of dict_simresults for epidemics only
# dict_epiincid[(code, sim number, T/C/A/D/S)] = cases over duration of epidemic

### other variables
# episize = number of nodes infected that defines an epidemic on the network (often defined as 5% of network size)

#####################

### import packages/modules ###
import random as rnd
import numpy as np
from collections import defaultdict
import zipfile
import bisect
from time import clock

# 7/15/14: functions rearranged in alphabetic order

####################################################
def atrisk_susceptibles(G, states):
	''' Return the list of susceptible nodes with at least one infected neighbor.
	'''
	return [u for u in states if states[u] == 's' and infected_neighbors(G, u, states) > 0]

####################################################
def calc_OR_from_list(dict_node_age, infected_nodelist):
	""" Calculate OR from list of infected nodes. Return OR as 1.0 if numerator or denominator are 0. 
	"""
	main(calc_OR_from_list)

	infected_child = sum([1 for node in infected_nodelist if dict_node_age[node] == '3'])
	infected_adult = sum([1 for node in infected_nodelist if dict_node_age[node] == '4'])
	
	# assign size of child and adult populations in network
	num_child, num_adult = child_adult_size(dict_node_age)
	
	# calculate incidence in network
	preval_child = infected_child/float(num_child)
	preval_adult = infected_adult/float(num_adult)
	
	# calculate OR if incid greater than 0 in both groups
	if (preval_child > 0 and preval_child < 1 and preval_adult > 0 and preval_adult < 1):
		OR = (preval_child/(1 - preval_child)) / (preval_adult/(1 - preval_adult))
	else:
		# not a number if 0 or division by 0
		OR = float('NaN')

	# return value
	return OR


####################################################
def calculate_beta(T, gamma, susc_node, dict_node_age, dict_age_susceptibility):
	""" Calculate beta for a given susceptible node based on T and a susceptibility value for the susceptible node. Susceptibility values range from 0 (not susceptible at all) to 1 (completely susceptible). Susceptibility values are assigned by age class. This calculation operates under the assumption that T = infectivity * susceptibility.
	"""
	main(calculate_beta)

	T_mod = T * dict_age_susceptibility[dict_node_age[susc_node]]
	beta_mod = (-T_mod * gamma)/(T_mod - 1)
	
	# return beta value based on modified T value
	return beta_mod

####################################################
def child_adult_size(dict_node_age):
	''' Return the number of child and adult nodes in the network. 
	'''
	main(child_adult_size)
	child_size = len([node for node in dict_node_age if dict_node_age[node] =='3'])
	adult_size = len([node for node in dict_node_age if dict_node_age[node] =='4'])
	return float(child_size), float(adult_size)

####################################################
def define_epi_time(dict_epiincid, b_or_s, align_proportion):
	""" For time-based epidemic simulations, identify the tstep at which the simulation reaches align_proportion proportion of cumulative infections for the epidemic. The goal is to align the epidemic trajectories at the point where the tsteps are comparable for each simulation. Returns dict where key = beta or s and value = list of tstep at which simulations that reached epidemic sizes attained 5% cumulative infections. 
	"""
	main(define_epi_time)
	
	dummykeys = [k for k in dict_epiincid if k[0] == b_or_s and k[2] == 'T']
	# dict_dummyalign_tstep[b_or_s] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	dict_dummyalign_tstep = defaultdict(list)
	
	for k in dummykeys:
		cum_infections = np.cumsum(dict_epiincid[k])
		prop_cum_infections = [float(item)/cum_infections[-1] for item in cum_infections]
		dict_dummyalign_tstep[b_or_s].append(bisect.bisect(prop_cum_infections, align_proportion))
	
	# change this value if you want the plots to align at a different tstep
	match_tstep = int(round(np.mean(dict_dummyalign_tstep[b_or_s]))) 
	
	return dict_dummyalign_tstep, match_tstep, dummykeys

###################################################
def epidemicsonly(dict_simresults, episize):
	''' Separate epidemics from simulations with an epidemic size less than or equal to episize.
	'''
	main(epidemicsonly)
	return dict((k, dict_simresults[k]) for k in dict_simresults if dict_simresults[k][2] > episize)

####################################################
def episim_age_time(G, dict_node_age, beta, gamma):
	''' Time-based simulation function for age-structured networks that returns the number of total infected individuals, infection time step per node, and recovery time step per node during the simulation. A graph, node-age class dictionary, and beta and gamma values must be provided.
	'''	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
		
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# number of nodes in graph
	N = int(G.order())
	
	# record infection timestep for patient zero
	# node 1 will be in column index 0 (nodes number from 1 to 10304)
	I_tstep_savelist = [0] * N
	I_tstep_savelist[int(p_zero)-1] = tstep 
	
	# create list to record recovery timesteps
	R_tstep_savelist = [0] * N

### simulation ###
	while infected_tstep:
		tstep += 1
				
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's']
		new_infected = [u for u in suscep_tstep if rnd.random() < (1- np.exp(-beta*infected_neighbors(G, u, states)))]
		for inf in new_infected:
			I_tstep_savelist[int(inf)-1] = tstep

		# I to R
		new_recovered = [v for v in infected_tstep if rnd.random() < gamma]
		for rec in new_recovered:
			R_tstep_savelist[int(rec) - 1] = tstep

		# update list of currently infected nodes for next tstep
		infected_tstep = infected_tstep + new_infected
		infected_tstep = [inf for inf in infected_tstep if inf not in new_recovered]
		# update states dictionary
		for new_i in new_infected:
			states[new_i] = 'i'
		for new_r in new_recovered:
			states[new_r] = 'r'
	
	# report total epidemic size
	recovered = [node for node in states if states[node] == 'r']

	### return data structures ###
	return len(recovered), I_tstep_savelist, R_tstep_savelist

####################################################
def episim_age_time_imm(G, dict_node_age, dict_immunity_mag, beta, gamma):
	''' Time-based simulation function for age-structured networks that returns the number of total infected individuals, infection time step per node, and recovery time step per node during the simulation. A graph, node-age class dictionary, node-pre-existing immunity dictionary, beta and gamma values must be provided. This simulation will be used generally to model the effects of heterogeneity in pre-existing immunity for the adult age class.
	'''	
	# initiate states dict
	states = dict([(node, 's') for node in G.nodes()]) # states[nodenumber] = state

	# time steps begin at 0
	tstep = 0
	# number of nodes in graph
	N = int(G.order())

	# node 1 will be in column index 0 (nodes number from 1 to 10304)
	I_tstep_savelist = [0] * N
	# create list to record recovery timesteps
	R_tstep_savelist = [0] * N

	# Choose seed infections
	p_zero, states, I_tstep_savelist = set_initial_seeds(G, states, I_tstep_savelist, 1)
	# list of currently infected nodes
	infected_tstep = p_zero[:]
	
### simulation ###
	while infected_tstep:
		tstep += 1
				
		# S to I
		# iterate through only susceptibles with infected neighbors
		suscep_tstep = atrisk_susceptibles(G, states)
		# list of new infections this time step
		new_infected = [u for u in suscep_tstep if rnd.random() < (1-dict_immunity_mag[u])*(1- np.exp(-beta*infected_neighbors(G, u, states)))]
		# update tracking dictionaries
		for inf in new_infected:
			states[inf] = 'i'
			I_tstep_savelist[int(inf)-1] = tstep

		# I to R
		# iterate through all infecteds that were infected prior to this time step
		new_recovered = [v for v in infected_tstep if rnd.random() < gamma]
		# update tracking dictionaries
		for rec in new_recovered:
			R_tstep_savelist[int(rec) - 1] = tstep
			states[rec] = 'r'

		# update list of currently infected nodes for next tstep
		infected_tstep = infected_tstep + new_infected
		infected_tstep = [inf for inf in infected_tstep if inf not in new_recovered]
	
	# report total epidemic size
	recovered = [node for node in states if states[node] == 'r']

	### return data structures ###
	return len(recovered), I_tstep_savelist, R_tstep_savelist


####################################################
def episim_age_time_realistic(G, dict_node_age, dict_age_params):
	''' Realistic time-based age-structured simulation function. Sigma, beta, and gamma may vary for children and adults at the same time.
	'''
	main(episim_age_time_realistic)
	
	# pull parameters from dict
	sigma_c, beta_c, gamma_c = dict_age_params['3']
	sigma_nc, beta_nc, gamma_nc = dict_age_params['4']
	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
	
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# number of nodes in graph
	N = int(G.order())
	
	# record infection timestep for patient zero
	# node 1 will be in column index 0 (nodes number from 1 to 10304)
	I_tstep_savelist = [0] * N
	I_tstep_savelist[int(p_zero)-1] = tstep 
	
	# create list to record recovery timesteps
	R_tstep_savelist = [0] * N

### simulation ###
	while infected_tstep:
		tstep += 1
		
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's']
		# dict_age_params[dict_node_age[u]][0] equals sigma of u
		
		# number of infected child neighbors for each susceptible 
		icn = [infected_child_neighbors(G, u, states, dict_node_age) for u in suscep_tstep]
		# number of infected non-child neighbors for each susceptible
		inn = [infected_nonchild_neighbors(G, u, states, dict_node_age) for u in suscep_tstep]
		# d_infnei[suscep node] = (# inf c_neighbors, # inf nc_neighbors)
		d_infnei = dict(zip(suscep_tstep, zip(icn, inn))) 
		# exclude values in dictionary where there are zero infected neighbors
		d_infnei_sub = dict((k, d_infnei[k]) for k in d_infnei if sum(d_infnei[k]) > 0)
		suscep_tstep_sub = [k for k in d_infnei_sub]
				
		new_infected = [u for u in suscep_tstep_sub if rnd.random() < dict_age_params[dict_node_age[u]][0] * (1- np.exp(-beta_c * d_infnei_sub[u][0]) * np.exp(-beta_nc * d_infnei_sub[u][1]))]

		for inf in new_infected:
			I_tstep_savelist[int(inf)-1] = tstep
		
		# I to R
		new_recovered_children = [v for v in infected_tstep if dict_node_age[v] == '3' and rnd.random() < gamma_c]
		new_recovered_nonchildren = [v for v in infected_tstep if dict_node_age[v] != '3' and rnd.random() < gamma_nc]
		new_recovered = new_recovered_children + new_recovered_nonchildren
		for rec in new_recovered:
			R_tstep_savelist[int(rec)-1] = tstep
		
		# update list of currently infected nodes for next tstep
		infected_tstep = infected_tstep + new_infected
		infected_tstep = [inf for inf in infected_tstep if inf not in new_recovered]
		# update states dictionary
		for new_i in new_infected:
			states[new_i] = 'i'
		for new_r in new_recovered:
			states[new_r] = 'r'
			
	# report total epidemic size
	recovered = [node for node in states if states[node] == 'r']
	
	### return data structures
	return len(recovered), I_tstep_savelist, R_tstep_savelist

####################################################
def episim_age_time_rec(G, dict_node_age, dict_age_recovery):
	''' Time-based age-structured simulation function that returns the number of total infected individuals, infection time step per node, and recovery time step per node during the simulation. Different recovery rates (1/infectious period) may be set for each age group in the main code and imported as pre-calculated betas in dict_age_recovery.
	'''
	main(episim_age_time_rec)
	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
	
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# number of nodes in graph
	N = int(G.order())
	
	# record infection timestep for patient zero
	# node 1 will be in column index 0 (nodes number from 1 to 10304)
	I_tstep_savelist = [0] * N
	I_tstep_savelist[int(p_zero)-1] = tstep 
	
	# create list to record recovery timesteps
	R_tstep_savelist = [0] * N

### simulation ###
	while infected_tstep:
		tstep += 1
				
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's']

		# number of infected non-child neighbors for each susceptible
		inei = [infected_neighbors(G, u, states) for u in suscep_tstep]
		# d_infnei[suscep node] = # infected neighbors
		d_infnei = dict(zip(suscep_tstep, inei)) 
		# exclude values in dictionary where there are zero infected neighbors
		d_infnei_sub = dict((k, d_infnei[k]) for k in d_infnei if d_infnei[k] > 0)
		suscep_tstep_sub = [k for k in d_infnei_sub]
				
		new_infected = [u for u in suscep_tstep_sub if rnd.random() < (1- np.exp(-dict_age_recovery[dict_node_age[u]][0] * d_infnei_sub[u]))]

		
		for inf in new_infected:
			I_tstep_savelist[int(inf)-1] = tstep

		# I to R
		new_recovered = [v for v in infected_tstep if rnd.random() < dict_age_recovery[dict_node_age[v]][1]]
		for rec in new_recovered:
			R_tstep_savelist[int(rec)-1] = tstep
		
		# update list of currently infected nodes for next tstep
		infected_tstep = infected_tstep + new_infected
		infected_tstep = [inf for inf in infected_tstep if inf not in new_recovered]
		# update states dictionary
		for new_i in new_infected:
			states[new_i] = 'i'
		for new_r in new_recovered:
			states[new_r] = 'r'
			
	# report total epidemic size
	recovered = [node for node in states if states[node] == 'r']

	### return data structures ###
	return len(recovered), I_tstep_savelist, R_tstep_savelist

####################################################
def episim_age_time_sickbehav(G, dict_node_age, beta, gamma, delay, cut, pop):
	''' Time-based simulation function for age-structured networks that returns the number of total infected individuals, infection time step per node, and recovery time step per node during the simulation. A graph, node-age class dictionary, and beta and gamma values must be provided. Infected children, adults, or both may lose 'cut' proportion of school or work contacts during infectious period after 'delay' days.
	'''	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()]) # states[nodenumber] = state
	states_tstep = dict([(node, 0) for node in G.nodes()]) # states_tstep[nodenumber] = tstep at which current state was achieved
	
	# # Randomly choose one node as patient zero
	# p_zero = rnd.choice(G.nodes()) 
	# states[p_zero] = 'i'
	# states_tstep[p_zero] = 0
	# infected_tstep = [p_zero]

	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	
	# number of nodes in graph
	N = int(G.order())

	# record infection timestep for patient zero
	# node 1 will be in column index 0 (nodes number from 1 to 10304)
	I_tstep_savelist = [0] * N
	# create list to record recovery timesteps
	R_tstep_savelist = [0] * N

	# Randomly choose seed infections
	p_zero = list(rnd.sample(G.nodes(), 1))
	for p in p_zero:
		states[p] = 'i'
		states_tstep[p] = tstep
		I_tstep_savelist[int(p)-1] = tstep
	
	infected_tstep = p_zero[:]
	

### simulation ###
	while infected_tstep:
		tstep += 1
				
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's' and infected_neighbors(G, u, states) > 0]

		new_infected = [u for u in suscep_tstep if rnd.random() < (1- np.exp(-beta*infected_neighbors_sickbehav(G, u, states, states_tstep, dict_node_age, delay, pop, cut, tstep)))]

		for inf in new_infected:
			I_tstep_savelist[int(inf)-1] = tstep

		# I to R
		new_recovered = [v for v in infected_tstep if rnd.random() < gamma]
		for rec in new_recovered:
			R_tstep_savelist[int(rec) - 1] = tstep

		# update list of currently infected nodes for next tstep
		infected_tstep = infected_tstep + new_infected
		infected_tstep = [inf for inf in infected_tstep if inf not in new_recovered]
		# update states and states_tstep dictionaries
		for new_i in new_infected:
			states[new_i] = 'i'
			states_tstep[new_i] = tstep
		for new_r in new_recovered:
			states[new_r] = 'r'
			states_tstep[new_r] = tstep
	
	# report total epidemic size
	recovered = [node for node in states if states[node] == 'r']

	### return data structures ###
	return len(recovered), I_tstep_savelist, R_tstep_savelist

	
####################################################
def episim_age_time_T(G, dict_node_age, beta, gamma, dict_age_betamodified):
	''' Time-based age-structured simulation function that returns the number of total infected individuals, infection time step per node, and recovery time step per node during the simulation. Different T multipliers may be set for each age group.
	'''	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
	
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# number of nodes in graph
	N = int(G.order())
	
	# record infection timestep for patient zero
	# node 1 will be in column index 0 (nodes number from 1 to 10304)
	I_tstep_savelist = [0] * N
	I_tstep_savelist[int(p_zero)-1] = tstep 
	
	# create list to record recovery timesteps
	R_tstep_savelist = [0] * N
	
	# beta vals
	b3, b4 = dict_age_betamodified['3'], dict_age_betamodified['4']

### simulation ###
	while infected_tstep:
		tstep += 1
				
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's']

		# beta is different for infected child and non-child neighbors

		# number of infected child neighbors for each susceptible 
		icn = [infected_child_neighbors(G, u, states, dict_node_age) for u in suscep_tstep]
		# number of infected non-child neighbors for each susceptible
		inn = [infected_nonchild_neighbors(G, u, states, dict_node_age) for u in suscep_tstep]
		# d_infnei[suscep node] = (# inf c_neighbors, # inf nc_neighbors)
		d_infnei = dict(zip(suscep_tstep, zip(icn, inn))) 
		# exclude values in dictionary where there are zero infected neighbors
		d_infnei_sub = dict((k, d_infnei[k]) for k in d_infnei if sum(d_infnei[k]) > 0)
		suscep_tstep_sub = [k for k in d_infnei_sub]
				
		new_infected = [u for u in suscep_tstep_sub if rnd.random() < (1- np.exp(-b3 * d_infnei_sub[u][0]) * np.exp(-b4 * d_infnei_sub[u][1]))]
		
		for inf in new_infected:
			I_tstep_savelist[int(inf)-1] = tstep

		# I to R
		new_recovered = [v for v in infected_tstep if rnd.random() < gamma]
		for rec in new_recovered:
			R_tstep_savelist[int(rec)-1] = tstep
		
		# update list of currently infected nodes for next tstep
		infected_tstep = infected_tstep + new_infected
		infected_tstep = [inf for inf in infected_tstep if inf not in new_recovered]
		# update states dictionary
		for new_i in new_infected:
			states[new_i] = 'i'
		for new_r in new_recovered:
			states[new_r] = 'r'
			
	# report total epidemic size
	recovered = [node for node in states if states[node] == 'r']

	### return data structures ###
	return len(recovered), I_tstep_savelist, R_tstep_savelist

####################################################
def episim_age_time_susc(G, dict_node_age, beta, gamma, dict_age_susceptibility):
	''' Time-based age-structured simulation function that returns the number of total infected individuals, infection time step per node, and recovery time step per node during the simulation. Different susceptibility values may be set for each age group.
	'''	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
	
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# number of nodes in graph
	N = int(G.order())
	
	# record infection timestep for patient zero
	# node 1 will be in column index 0 (nodes number from 1 to 10304)
	I_tstep_savelist = [0] * N
	I_tstep_savelist[int(p_zero)-1] = tstep 
	
	# create list to record recovery timesteps
	R_tstep_savelist = [0] * N

### simulation ###
	while infected_tstep:
		tstep += 1
				
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's']
		
		# number of infected non-child neighbors for each susceptible
		inei = [infected_neighbors(G, u, states) for u in suscep_tstep]
		# d_infnei[suscep node] = # infected neighbors
		d_infnei = dict(zip(suscep_tstep, inei)) 
		# exclude values in dictionary where there are zero infected neighbors
		d_infnei_sub = dict((k, d_infnei[k]) for k in d_infnei if d_infnei[k] > 0)
		suscep_tstep_sub = [k for k in d_infnei_sub]
				
		new_infected = [u for u in suscep_tstep_sub if rnd.random() < dict_age_susceptibility[dict_node_age[u]] * (1- np.exp(-beta * d_infnei_sub[u]))] # 3/22/14 sigma was in the wrong place


		for inf in new_infected:
			I_tstep_savelist[int(inf)-1] = tstep

		# I to R
		new_recovered = [v for v in infected_tstep if rnd.random() < gamma]
		for rec in new_recovered:
			R_tstep_savelist[int(rec)-1] = tstep
		
		# update list of currently infected nodes for next tstep
		infected_tstep = infected_tstep + new_infected
		infected_tstep = [inf for inf in infected_tstep if inf not in new_recovered]
		# update states dictionary
		for new_i in new_infected:
			states[new_i] = 'i'
		for new_r in new_recovered:
			states[new_r] = 'r'
			
	# report total epidemic size
	recovered = [node for node in states if states[node] == 'r']

	### return data structures ###
	return len(recovered), I_tstep_savelist, R_tstep_savelist

####################################################
def filter_time_points2(tot_incidlist, incl_min, incl_max):
	""" pairs with reg time sim only (11/6/13) In a time-based epidemic simulation, filter OR data points that may not be valid due to the small number of infected individuals at the beginning and end of an epidemic. For each simulation, define this period as the time period between which the percentage of infecteds is incl_min (float between 0 and 1) and incl_max (float between 0 and 1) percent of the cumulative epidemic size for that simulation. 
	"""

	# calculate cumulative number of cases at each time step
	cum_ct = np.cumsum(tot_incidlist)
	# generate list of time steps where cumulative number of cases is between incl_min and incl_max proportion of total cases in epidemic
	incl_tsteps = [i for i, ct in enumerate(cum_ct) if ct/float(max(cum_ct)) > incl_min and ct/float(max(cum_ct)) < incl_max]
	
	if incl_tsteps:
		return min(incl_tsteps), max(incl_tsteps)
	else:
		return 0, 0

####################################################
def immune_nodes_proportion(G, dict_node_age, proportion):
	''' Choose a random proportion of the adult population to have pre-existing immunity.
	'''
	rnd.seed(107)
	_, num_adults = child_adult_size(dict_node_age)
	number_immune = int(round(num_adults * proportion, 0))
	adult_nodes = [node for node in dict_node_age if dict_node_age[node] == '4']
	print "number of immune adults", number_immune
	immune_nodes = rnd.sample(adult_nodes, number_immune)

	# perform check that immune_nodes returned only adults
	print "all immune_nodes are adults: ", len(immune_nodes) == sum([1 for i in immune_nodes if dict_node_age[i] == '4'])

	return immune_nodes

####################################################
def immune_nodes_simulation(G, dict_node_age, beta, gamma):
	''' Use results of an epidemic simulation with the same transmission parameters in order to choose adult nodes with pre-existing immunity in the subsequent simulations.
	'''
	rnd.seed(107)
	# create binary lists to indicate children and adults
	adults = [1 if dict_node_age[str(node)] == '4' else 0 for node in xrange(1, int(N) + 1)]
	# perform a single simulation to choose adults with immunity
	episize, I_tstep_savelist, R_tstep_savelist = episim_age_time(G, dict_node_age, beta, gamma)
	# binary list of adults infected in simulation
	infected_adults = [1 if a*t > 0 else 0 for a, t in zip(adults, R_tstep_savelist)]
	# list of all node IDs in network
	indexes = range(1, len(infected_adults)+1)
	# list of infected adult node IDs
	immune_nodes = [a*i for a,i in zip(infected_adults, indexes) if a]

	# perform check that immune_nodes returned only adults
	print "all immune_nodes are adults: ", len(immune_nodes) == sum([1 for i in immune_nodes if dict_node_age[i] == '4'])

	return immune_nodes

####################################################
def immunity_magnitude_single(G, immune_nodes, immune_value):
	''' Assign magnitude of pre-existing immunity 0 to 1 for each node. A single immune_value magnitude of pre-existing immunity exists, and it represents the average level of immunity in the subpopulation defined by the list immune_nodes.
	'''
	dict_immunity_mag = dict(zip(G.nodes(), [0]*G.order()))
	for node in immune_nodes:
		dict_immunity_mag[node] = immune_value
	
	return dict_immunity_mag

####################################################
def immunity_magnitude_uniform(G, immune_nodes, max_immune_value):
	''' Assign magnitude of pre-existing immunity 0 to 1 for each node. Draw immunity values from a uniform distribution with the maximum immunity level set at max_immune_value.
	'''
	dict_immunity_mag = dict(zip(G.nodes(), [0]*G.order()))
	rnd.seed(88)
	for node in immune_nodes:
		dict_immunity_mag[node] = rnd.uniform(0,max_immune_value)
	
	return dict_immunity_mag

####################################################
def infected_child_neighbors(G, node, states, dict_node_age):
	""" Calculate the number of infected child neighbors for the node. 
	"""
	return sum([1 for node_i in G.neighbors(node) if states[node_i] == 'i' and dict_node_age[node_i] == '3'])

####################################################
def infected_nonchild_neighbors(G, node, states, dict_node_age):
	""" Calculate the number of infected non-child neighbors for the node. 
	"""
	return sum([1 for node_i in G.neighbors(node) if states[node_i] == 'i' and dict_node_age[node_i] != '3'])

####################################################
def infected_neighbors(G, node, states):
	""" Calculate the number of infected neighbors for the node. 
	"""
	return sum([1 for node_i in G.neighbors(node) if states[node_i] == 'i']) 

####################################################
def infected_neighbors_sickbehav(G, node, states, states_tstep, dict_node_age, delay, pop, cut, tstep):
	""" Calculate the number of infected neighbors for the susceptible node. Infected neighbors cannot infect work or school contacts after 'delay' number of days if they are adults or children, respectively. 'pop' parameter designates whether the simulation cuts these contacts for infected adults, children, or both. 'cut' parameter designates the proportion of contacts that are cut after 'delay' number of days has passed for infected adults/children with work/school contacts.
	"""
	pop_options = ['A', 'C', 'AC', 'CA']
	active_nonschoolwork, active_nonschool, active_schoolnonchild, active_schoolchilddelay, active_schoolchildpostdelay, active_nonwork, active_worknonadult, active_workadultdelay, active_workadultpostdelay = [],[],[],[],[],[],[],[],[]
	active_neighbors = [nei for nei in G.neighbors(node) if states[nei] == 'i']
	
	if pop not in pop_options:
		print 'pop parameter error'

	elif pop == 'C':
		active_nonschool = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] != 'S']
		active_schoolnonchild = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'S' and dict_node_age[nei] != '3'] # school contact where infected neighbor is not a child
		active_schoolchilddelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'S' and dict_node_age[nei] == '3' and tstep - states_tstep[nei] < delay] # schoolchild may always infect at school before "delay" number of days has passed during infectious period 
		active_schoolchildpostdelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'S' and dict_node_age[nei] == '3' and tstep - states_tstep[nei] >= delay and rnd.random() > cut] # schoolchild may infect at school after "delay" number of days "1-cut" proportion of the time; this is an approximation because the identity of the individuals that are staying home is not fixed

		if sorted(active_nonschool + active_schoolnonchild + active_schoolchilddelay + active_schoolchildpostdelay) != sorted(set(active_nonschool + active_schoolnonchild + active_schoolchilddelay + active_schoolchildpostdelay)):
			print "pop C error" # check that neighbors were not double counted
	
	elif pop == 'A':
		active_nonwork = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] != 'W']
		active_worknonadult = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'W' and dict_node_age[nei] != '4']
		active_workadultdelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'W' and dict_node_age[nei] == '4' and tstep - states_tstep[nei] < delay]
		active_workadultpostdelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'W' and dict_node_age[nei] == '4' and tstep - states_tstep[nei] >= delay and rnd.random() > cut] 

		if sorted(active_nonwork + active_worknonadult + active_workadultdelay + active_workadultpostdelay) != sorted(set(active_nonwork + active_worknonadult + active_workadultdelay + active_workadultpostdelay)):
			print "pop A error"

	else: # if pop == 'AC' or 'CA'
		active_nonschoolwork = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] not in list('WS')] # neither school nor work infected contacts because conditions apply to both of those edge types
		active_schoolnonchild = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'S' and dict_node_age[nei] != '3'] 
		active_schoolchilddelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'S' and dict_node_age[nei] == '3' and tstep - states_tstep[nei] < delay] 
		active_schoolchildpostdelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'S' and dict_node_age[nei] == '3' and tstep - states_tstep[nei] >= delay and rnd.random() > cut]
		active_worknonadult = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'W' and dict_node_age[nei] != '4']
		active_workadultdelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'W' and dict_node_age[nei] == '4' and tstep - states_tstep[nei] < delay]
		active_workadultpostdelay = [nei for nei in active_neighbors if G.get_edge_data(node, nei)['place'] == 'W' and dict_node_age[nei] == '4' and tstep - states_tstep[nei] >= delay and rnd.random() > cut] 

		# print 
		# print 'node, node-age, nei, nei-age', node, dict_node_age[node],active_neighbors, [dict_node_age[nei] for nei in active_neighbors]
		# print [G.get_edge_data(node, nei)['place'] for nei in active_neighbors]
		# print 'nonschoolwork', active_nonschoolwork
		# print 'schoolnonchild', active_schoolnonchild
		# print 'schoolchilddelay', active_schoolchilddelay
		# print 'schoolchildpostdelay', active_schoolchildpostdelay
		# print 'worknonadult', active_worknonadult
		# print 'workadultdelay', active_workadultdelay
		# print 'workadultpostdelay', active_workadultpostdelay

		if sorted(active_nonschoolwork + active_worknonadult + active_workadultdelay + active_workadultpostdelay + active_schoolnonchild + active_schoolchilddelay + active_schoolchildpostdelay) != sorted(set(active_nonschoolwork + active_worknonadult + active_workadultdelay + active_workadultpostdelay + active_schoolnonchild + active_schoolchilddelay + active_schoolchildpostdelay)):
			print "pop AC error"

	neighbor_ct = len(active_nonschoolwork) + len(active_nonschool) + len(active_schoolnonchild) + len(active_schoolchilddelay) + len(active_schoolchildpostdelay) + len(active_nonwork) + len(active_worknonadult) + len(active_workadultdelay) + len(active_workadultpostdelay)	

	return float(neighbor_ct)

####################################################
def mean_retro_OR(dict_epiincid, dict_epiOR, code, min_percent, max_percent):
	''' Return a list of the mean retro OR values for all simulations with a single set of parameters. The retrospective period is defined as the time steps that correspond with the 'min_percent' and 'max_percent' cumulative incidence percentages. These percentages were derived from the average cumulative incidence percentage of the retrospective period in the SDI data across all 8 seasons in the study.
	'''
	main(mean_retro_OR)

	keys = [k for k in dict_epiincid if k[0] == code and k[2] == 'T']
	mn_retro_OR_ls = []
	for k in keys:
		tmin, tmax = filter_time_points2(dict_epiincid[k], min_percent, max_percent)
		dummyarray = dict_epiOR[(k[0], k[1])][tmin:tmax+1]
		mn_retro_OR = np.mean(np.ma.masked_array(dummyarray, np.isnan(dummyarray))) # calculate mean without nan values
		mn_retro_OR_ls.append(mn_retro_OR)

	# report mean list without nan values
	mn_retro_OR_ls_masked = np.ma.masked_array(mn_retro_OR_ls, np.isnan(mn_retro_OR_ls))

	return mn_retro_OR_ls_masked

####################################################
def OR_sim(numbersims, dict_epiresults, Tmult, childsize, adultsize):
	""" Calculates overall OR for an entire simulation based on dict_epiresults.
	"""
	main(OR_sim)
	
	# dict_epiresults[(Tmult, simnumber)] = (episize, c_episize, a_episize)
	cAR = [dict_epiresults[key][1]/childsize for key in dict_epiresults if key[0] == Tmult]
	aAR = [dict_epiresults[key][2]/adultsize for key in dict_epiresults if key[0] == Tmult]
	# 3/28/14: there was an order of operations issue (parentheses) in the previous list comprehension, but this is now fixed
	OR_sim_ls = [(c/(1-c))/(a/(1-a)) for c, a in zip(cAR, aAR)]

	return OR_sim_ls

####################################################
def perc_age_gen(G, dict_node_age, T):
	''' General age-structured percolation function that returns the number of children, adults, and total infected individuals during the simulation.
	'''
	main(perc_age_gen)
	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	p_zero = rnd.choice(G.nodes()) # Randomly choose one node as patient zero
	states[p_zero] = 'i'
	infected    = [p_zero] 
	recovered   = [] 
	
	# simulation
	while infected:
		v = infected.pop(0)
		for u in G.neighbors(v):
			# make sure node u is still susceptible
			if states[u] == 's' and rnd.random() < T: 
				states[u] = 'i'
				infected.append(u)
		states[v] = 'r'
		recovered.append(v)
	
	# return binary list of ordered number nodes that were infected
	binary_list = [0 for n in xrange(G.order())]
	for node in recovered:
		binary_list[int(node) - 1] = 1
	
	# infected children
	rec_child_n = float(len([node for node in recovered if dict_node_age[node] == '3']))
	# infected adults
	rec_adult_n = float(len([node for node in recovered if dict_node_age[node] == '4']))
	
	# return OR for simulation
	ORval = calc_OR_from_list(dict_node_age, recovered)
	
	# number recovered in children and adults and total
	return rec_child_n, rec_adult_n, len(recovered), binary_list, ORval

####################################################
def perc_age_vaxeff(G, dict_node_age, T, vaxcov_scen, susceptibility):
	''' Age-structured percolation function with random vaccination with low (0.245) and high (0.455) coverage scenarios. Susceptibility (1 - vax efficacy) is set outside the function and vaccinated individuals have a probability of T * susceptibility of infection. The function returns the number of children, adults, total infected, and total vaccinated during th simulation. Susceptibility of 1 means fully susceptible.
	'''
	main(perc_age_vaxeff)

	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	p_zero = rnd.choice(G.nodes()) # Randomly choose one node as patient zero
	states[p_zero] = 'i'
	infected = [p_zero] 
	recovered = [] 
	vaccinated = []
	
	# define two vax coverage scenarios
	if vaxcov_scen == 'Low':
		vaxcov = 0.245
	elif vaxcov_scen == 'High':
		vaxcov = 0.455
	else:
		print "vaccine coverage scenario error"
	
	# identify vaccinated nodes
	for n in G.nodes():
		if rnd.random() < vaxcov and n != p_zero:
			states[n] = 'v'
			vaccinated.append(n)
	
	# simulation
	while infected:
		v = infected.pop(0)
		for u in G.neighbors(v):
			if states[u] == 's' and rnd.random() < T:
				states[u] = 'i'
				infected.append(u)
			elif states[u] == 'v' and rnd.random() < T * susceptibility:
				states[u] = 'i'
				infected.append(u)
		states[v] = 'r'
		recovered.append(v)
	
	# return binary list of ordered number nodes that were infected
	binary_list = [0 for n in np.arange(G.order())]
	for node in recovered:
		binary_list[int(node) - 1] = 1
	
	# infected children
	rec_child_n = float(len([node for node in recovered if dict_node_age[node] == '3']))
	# infected adults
	rec_adult_n = float(len([node for node in recovered if dict_node_age[node] == '4']) )
	
	# number recovered children, adults, and total pop, and number of initially vaccinated individuals
	return rec_child_n, rec_adult_n, len(recovered), len(vaccinated), binary_list


####################################################
def recreate_epidata(I_filename, R_filename, zipname, b_or_s, epi_size, child_nodes, adult_nodes, dict_epiincid, dict_epiOR, dict_epiresults, dict_epiAR, dict_epiOR_filt):
	""" For time-based epidemic simulations, recreate total and age-specific incidence, OR, and general sim results from simulation output file of timestep at which each node got infected, where column indexes are node IDs minus one and rows are simulation results. child_nodes and adult_nodes are binary lists indicating whether the nodeID is a child/adult or not. Function returns five dictionaries for epidemic simulations: incidence rate by group, OR, total sim results, attack rate (number of infections/pop size), and OR filtered to time steps where 5-95% of cumulative incidence takes place. This function may be used for time sims for T, suscep, and recovery designations.
	"""
	main(recreate_epidata)
	
	# filter time point parameters, 5% and 95% cumulative incidence
	incl_min, incl_max = 0.05, 0.95
	
	# initialize Itstep and Rtstep dictionaries
	# d_Itstep[(beta, simnumber)] = [infection tstep for node 1, infection tstep for node 2, ...]
	dict_Itstep = defaultdict(list)
	# d_Rtstep[(beta, simnumber)] = [recovery tstep for node 1, recovery tstep for node 2, ...]
	dict_Rtstep = defaultdict(list)
		
	# import output file from ziparchive (5 sec)
	with zipfile.ZipFile(zipname, 'r') as zf:
		Itstep = zf.open(I_filename)
		Rtstep = zf.open(R_filename)
		simnumber = 0
		for line1, line2 in zip(Itstep, Rtstep):
			infected_tsteps = line1.split(', ')
			recovered_tsteps = line2.split(', ')
			# insert epidemics only into dicts
			if sum([0 if not int(tstep) else 1 for tstep in infected_tsteps]) > epi_size:
				dict_Itstep[(b_or_s, simnumber)] = [int(t) for t in infected_tsteps]
				dict_Rtstep[(b_or_s, simnumber)] = [int(t) for t in recovered_tsteps]
			simnumber += 1
		
	# total population size
	N = float(len(child_nodes))
	# size of child and adult populations
	children, adults = float(sum(child_nodes)), float(sum(adult_nodes))

	for b_or_s, simnum in dict_Itstep:
		# initiate lists to store incidence in different age groups where list index is tstep
		# skip counting patient zero in incidence and attack rate for convenience
		tot_AR, c_AR, a_AR = [0],[0],[0]
		tot_incid, c_incid, a_incid = [0], [0], [0]
		# child infection tsteps only, all else nan
		c_tsteps = [c*t for c, t in zip(child_nodes, dict_Itstep[(b_or_s, simnum)])]
		# adult infection tsteps only, all else nan
		a_tsteps = [a*t for a, t in zip(adult_nodes, dict_Itstep[(b_or_s, simnum)])]
		
		# count number of new infections (0.15 sec)
		# tstep 0 is skipped
		for I_tstep in xrange(1, max(dict_Itstep[(b_or_s, simnum)]) + 1):
			# number of new infections per 100 in total popn at I_tstep
			dummy_tot_inf = len([inf for inf in dict_Itstep[(b_or_s, simnum)] if inf == I_tstep])
			# number of new child infections per 100 in child popn at I_tstep
			dummy_c_inf = len([inf for inf in c_tsteps if inf == I_tstep])
			# number of new adult infections per 100 in adult pop at I_tstep
			dummy_a_inf = len([inf for inf in a_tsteps if inf == I_tstep])
			
			# append attack rate of total, child, and adult infections in respective total, child, and adult populations to list where index represents tstep
			# use this to calculate OR
			tot_AR.append(dummy_tot_inf/N)
			c_AR.append(dummy_c_inf/children)
			a_AR.append(dummy_a_inf/adults)
			
			# append incidence of total, child, and adult infections as raw counts
			# use this for dict_epiresults
			tot_incid.append(dummy_tot_inf)
			c_incid.append(dummy_c_inf)
			a_incid.append(dummy_a_inf)
		
		# recover attack rates over time in dictionary for each sim
		dict_epiAR[(b_or_s, simnum, 'T')] = tot_AR
		dict_epiAR[(b_or_s, simnum, 'C')] = c_AR
		dict_epiAR[(b_or_s, simnum, 'A')] = a_AR
				
		# recover new cases over time in dictionary for each sim as raw counts
		dict_epiincid[(b_or_s, simnum, 'T')] = tot_incid
		dict_epiincid[(b_or_s, simnum, 'C')] = c_incid
		dict_epiincid[(b_or_s, simnum, 'A')] = a_incid
		
		# recover OR
		# if attack rate in children or adults is 0 or 1, float(nan) is returned because the OR cannot be appropriately calculated
		OR_list = [((c/(1-c))/(a/(1-a))) if c and c < 1 and a and a < 1 else float('nan') for c, a in zip(c_AR, a_AR)]
		dict_epiOR[(b_or_s, simnum)] = OR_list
		
		# recover general epi results
		dict_epiresults[(b_or_s, simnum)] = (sum(tot_incid), sum(c_incid), sum(a_incid))
		
		# filter time points where current outbreak size is too small (just beginning or towards the end of the epidemic)
		min_tstep, max_tstep = filter_time_points2(tot_incid, incl_min, incl_max)
		
		dict_epiOR_filt[(b_or_s, simnum)] = [float('nan') if (num < min_tstep or num > max_tstep) else OR for num, OR in enumerate(OR_list)]

	return dict_epiincid, dict_epiOR, dict_epiresults, dict_epiAR, dict_epiOR_filt

####################################################
def recreate_epidata2(I_filename, R_filename, zipname, b_or_s, epi_size, child_nodes, adult_nodes, toddler_nodes, senior_nodes, dict_epiincid, dict_epiOR, dict_epiresults, dict_epiAR, dict_epiOR_filt):
	""" For time-based epidemic simulations, recreate total and age-specific incidence, OR, and general sim results from simulation output file of timestep at which each node got infected, where column indexes are node IDs minus one and rows are simulation results. child_nodes, adult_nodes, toddler_nodes, and senior_nodes are binary lists indicating whether the nodeID is a child/adult/toddler/senior or not. Toddlers and seniors are considered high risk groups. Function returns five dictionaries for epidemic simulations: incidence rate by group, OR, total sim results, attack rate (number of infections/pop size), and OR filtered to time steps where 5-95% of cumulative incidence takes place. This function may be used for time sims for T, suscep, and recovery designations.
	"""
	main(recreate_epidata2)
	
	# filter time point parameters, 5% and 95% cumulative incidence
	incl_min, incl_max = 0.05, 0.95
	
	# initialize Itstep and Rtstep dictionaries
	# d_Itstep[(beta, simnumber)] = [infection tstep for node 1, infection tstep for node 2, ...]
	dict_Itstep = defaultdict(list)
	# d_Rtstep[(beta, simnumber)] = [recovery tstep for node 1, recovery tstep for node 2, ...]
	dict_Rtstep = defaultdict(list)
		
	# import output file from ziparchive (5 sec)
	with zipfile.ZipFile(zipname, 'r') as zf:
		Itstep = zf.open(I_filename)
		Rtstep = zf.open(R_filename)
		simnumber = 0
		for line1, line2 in zip(Itstep, Rtstep):
			infected_tsteps = line1.split(', ')
			recovered_tsteps = line2.split(', ')
			# insert epidemics only into dicts
			if sum([0 if not int(tstep) else 1 for tstep in infected_tsteps]) > epi_size:
				dict_Itstep[(b_or_s, simnumber)] = [int(t) for t in infected_tsteps]
				dict_Rtstep[(b_or_s, simnumber)] = [int(t) for t in recovered_tsteps]
			simnumber += 1
		
	# total population size
	N = float(len(child_nodes))
	# size of child and adult populations
	children, adults = float(sum(child_nodes)), float(sum(adult_nodes))

	for b_or_s, simnum in dict_Itstep:
		# initiate lists to store incidence in different age groups where list index is tstep
		# skip counting patient zero in incidence and attack rate for convenience
		tot_AR, c_AR, a_AR = [0],[0],[0]
		# total, children, adult, toddler, senior
		tot_incid, c_incid, a_incid, d_incid, s_incid = [0], [0], [0], [0], [0]
		# child infection tsteps only, all else nan
		c_tsteps = [c*t for c, t in zip(child_nodes, dict_Itstep[(b_or_s, simnum)])]
		# adult infection tsteps only, all else nan
		a_tsteps = [a*t for a, t in zip(adult_nodes, dict_Itstep[(b_or_s, simnum)])]
		# toddler infection tsteps only, all else nan
		d_tsteps = [d*t for d, t in zip(toddler_nodes, dict_Itstep[(b_or_s, simnum)])]
		# senior infection tsteps only, all else nan
		s_tsteps = [s*t for s, t in zip(senior_nodes, dict_Itstep[(b_or_s, simnum)])]
		
		# count number of new infections (0.15 sec)
		# tstep 0 is skipped
		for I_tstep in xrange(1, max(dict_Itstep[(b_or_s, simnum)]) + 1):
			# number of new infections in total popn at I_tstep
			dummy_tot_inf = len([inf for inf in dict_Itstep[(b_or_s, simnum)] if inf == I_tstep])
			# number of new child infections in child popn at I_tstep
			dummy_c_inf = len([inf for inf in c_tsteps if inf == I_tstep])
			# number of new adult infections in adult pop at I_tstep
			dummy_a_inf = len([inf for inf in a_tsteps if inf == I_tstep])
			# number of new toddler infections in toddler pop at I_tstep
			dummy_d_inf = len([inf for inf in d_tsteps if inf == I_tstep])
			# number of new toddler infections in senior pop at I_tstep
			dummy_s_inf = len([inf for inf in s_tsteps if inf == I_tstep])
			
			# append attack rate of total, child, and adult infections in respective total, child, and adult populations to list where index represents tstep
			# use this to calculate OR
			tot_AR.append(dummy_tot_inf/N)
			c_AR.append(dummy_c_inf/children)
			a_AR.append(dummy_a_inf/adults)
			
			# append incidence of total, child, and adult infections as raw counts
			tot_incid.append(dummy_tot_inf)
			c_incid.append(dummy_c_inf)
			a_incid.append(dummy_a_inf)
			d_incid.append(dummy_d_inf)
			s_incid.append(dummy_s_inf)
		
		# recover attack rates over time in dictionary for each sim
		dict_epiAR[(b_or_s, simnum, 'T')] = tot_AR
		dict_epiAR[(b_or_s, simnum, 'C')] = c_AR
		dict_epiAR[(b_or_s, simnum, 'A')] = a_AR
				
		# recover new cases over time in dictionary for each sim as raw counts
		dict_epiincid[(b_or_s, simnum, 'T')] = tot_incid
		dict_epiincid[(b_or_s, simnum, 'C')] = c_incid
		dict_epiincid[(b_or_s, simnum, 'A')] = a_incid
		dict_epiincid[(b_or_s, simnum, 'D')] = d_incid
		dict_epiincid[(b_or_s, simnum, 'S')] = s_incid

		
		# recover OR
		# if attack rate in children or adults is 0 or 1, float(nan) is returned because the OR cannot be appropriately calculated
		OR_list = [((c/(1-c))/(a/(1-a))) if c and c < 1 and a and a < 1 else float('nan') for c, a in zip(c_AR, a_AR)]
		dict_epiOR[(b_or_s, simnum)] = OR_list # round(b_or_s, 1) may be necessary for key[0]
		
		# recover general epi results
		dict_epiresults[(b_or_s, simnum)] = (sum(tot_incid), sum(c_incid), sum(a_incid))
		
		# filter time points where current outbreak size is too small (just beginning or towards the end of the epidemic)
		min_tstep, max_tstep = filter_time_points2(tot_incid, incl_min, incl_max)
		
		dict_epiOR_filt[(b_or_s, simnum)] = [float('nan') if (num < min_tstep or num > max_tstep) else OR for num, OR in enumerate(OR_list)]

	return dict_epiincid, dict_epiOR, dict_epiresults, dict_epiAR, dict_epiOR_filt

####################################################
def rnd_vax(G, vaxcov):
	''' Implement a random vax strategy on a graph where vax efficacy is 100% and vax coverage is vaxcov%.
	'''
	main(rnd_vax)
	for n in G.nodes():
		if rnd.random() < vaxcov:
			G.remove_edges_from(G.edges(n))

####################################################
def set_initial_seeds(G, states, I_tstep_savelist, numseeds):
	''' Set initial conditions for simulation. Choose numseeds initial infections as patient zero(es) and update tracking dictionaries to reflect these infections.
	'''
	p_zero = list(rnd.sample(G.nodes(), numseeds))
	for p in p_zero:
		states[p] = 'i'
		# record infection tstep for seed infections
		I_tstep_savelist[int(p)-1] = 0 

	return p_zero, states, I_tstep_savelist


##############################################
##############################################
# footer

def main(function):
	print 'Running', __name__, function.__name__

if __name__ == '__main__':
	print 'Executed from the command line'
	main()












