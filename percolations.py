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
# dict_epiincid[(T, 'C' or 'A')] = child or adult attack rate

### other variables
# episize = number of nodes infected that defines an epidemic on the network (often defined as 5% of network size)

#####################

### import packages/modules ###
import random as rnd
import numpy as np
from collections import defaultdict
import zipfile
import bisect
import math



# functions
####################################################
# number of child and adult nodes in network
def child_adult_size(dict_node_age):
	''' Return the number of child and adult nodes in the network. 
	'''
	child_size = len([node for node in dict_node_age if dict_node_age[node] =='3'])
	adult_size = len([node for node in dict_node_age if dict_node_age[node] =='4'])
	return float(child_size), float(adult_size)

###################################################
# separate epidemics from all results
def epidemicsonly(dict_simresults, episize):
	''' From a dictionary of all results, return a dictionary subset of the epidemic results only. 
	'''
	return dict((k, dict_simresults[k]) for k in dict_simresults if dict_simresults[k][2] > episize)

####################################################
# random vax strategy with vaxcov% coverage and 100% efficacy
def rnd_vax(G, vaxcov):
	''' Implement a random vax strategy on a graph where vax efficacy is 100% and vax coverage is vaxcov%.
	'''
	for n in G.nodes():
		if rnd.random() < vaxcov:
			G.remove_edges_from(G.edges(n))

####################################################
# age structured percolation function with random vax strategy with vaxcov% coverage and defined susceptibility multiplier, where 1 is fully susceptible
def perc_age_vaxeff(G, dict_node_age, T, vaxcov_scen, susceptibility):
	''' Age-structured percolation function with random vaccination with low (0.245) and high (0.455) coverage scenarios. Susceptibility (1 - vax efficacy) is set outside the function and vaccinated individuals have a probability of T * susceptibility of infection. The function returns the number of children, adults, total infected, and total vaccinated during th simulation.
	'''

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
# general age-structured percolation function 
def perc_age_gen(G, dict_node_age, T):
	''' General age-structured percolation function that returns the number of children, adults, and total infected individuals during the simulation.
	'''
	
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
# time-based age-structured simulation function 
def episim_age_time(G, dict_node_age, beta, gamma):
	''' Time-based age-structured simulation function that returns the number of children, adults, total infected individuals, child to adult attack rate odds ratio at each time step, list of total incidence at each time step, list of total prevalence at each timestep, odds ratio at each time step with a substantial number of infections, OR for entire simulation, infection time step per node, and recovery time step per node during the simulation. 
	'''
	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
	
	# define cumulative percentiles that filter OR time points
	incl_min, incl_max = 0.1, 0.9
	
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# record infection timestep for patient zero
	I_tstep_savelist = [float('nan') for n in G.nodes()]
	I_tstep_savelist[int(p_zero)-1] = tstep # node 1 will be in column index 0 (nodes number from 1 to 10304)
	
	# create list to record recovery timesteps
	R_tstep_savelist = [float('nan') for n in G.nodes()]

	# tot_incidlist is a list of total incidence for each time step
	tot_incidlist = [len(infected_tstep)]
	# tot_prevallist is a list of total prevalence for each time step
	tot_prevallist = [len(infected_tstep)]
	# ORlist is a list of ORs for each time step
	ORlist = []

### simulation ###
	while infected_tstep:
		tstep += 1
				
		# count to keep track of incidence
		incid_ct = 0
				
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's']
		new_infected = [u for u in suscep_tstep if rnd.random() < (1- np.exp(-beta*infected_neighbors(G, u, states)))]
		for inf in new_infected:
			I_tstep_savelist[int(u)-1] = tstep

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
		
		incid_ct = len(new_infected)
		
### metrics by time step ###
		# 1) track total incidence for each time step
		tot_incidlist.append(incid_ct)
		
		# 2) track total prevalence for each time step
		tot_prevallist.append(len(infected_tstep))
			
		# 3) return a list of ORs for each time step
		OR = calc_OR_from_list(dict_node_age, infected_tstep)
		ORlist.append(OR)

### metrics over entire simulation ###
	# 1) For list of nodes that were infected over the entire simulation, look at indexes in I_tstep_savelist or R_tstep_savelist that are not float('nan')s. Remember that index = node number - 1
	recovered = [node for node in states if states[node] == 'r']
	# 2) infected children
	rec_child_n = float(len([node for node in recovered if dict_node_age[node] == '3']))
	# 3) infected adults
	rec_adult_n = float(len([node for node in recovered if dict_node_age[node] == '4']))

	# 4) OR dict filtered to include only time points with a substantial number of infections
	min_tstep, max_tstep = filter_time_points(tot_incidlist, recovered, incl_min, incl_max)
	filtered_ORlist = [float('NaN') if (num < min_tstep or num > max_tstep) else OR for num, OR in enumerate(ORlist)]

	# 5) total OR value
	ORval_total = calc_OR_from_list(dict_node_age, recovered)


	### return data structures ###
	return rec_child_n, rec_adult_n, len(recovered), ORlist, tot_incidlist, tot_prevallist, filtered_ORlist, ORval_total, I_tstep_savelist, R_tstep_savelist


####################################################
# time-based age-structured simulation function with varying susceptibilities by age group
def episim_age_time_susc(G, dict_node_age, T, gamma, dict_age_susceptibility):
	''' Time-based age-structured simulation function that returns the number of children, adults, total infected individuals, child to adult attack rate odds ratio at each time step, list of total incidence at each time step, list of total prevalence at each timestep, odds ratio at each time step with a substantial number of infections, OR for entire simulation, infection time step per node, and recovery time step per node during the simulation. Different susceptibility values may be set for each age group.
	'''
	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
	
	# define cumulative percentiles that filter OR time points
	incl_min, incl_max = 0.1, 0.9
	
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# record infection timestep for patient zero
	I_tstep_savelist = [float('nan') for n in G.nodes()]
	I_tstep_savelist[int(p_zero)-1] = tstep # node 1 will be in column index 0 (nodes number from 1 to 10304)
	
	# create list to record recovery timesteps
	R_tstep_savelist = [float('nan') for n in G.nodes()]

	# tot_incidlist is a list of total incidence for each time step
	tot_incidlist = [len(infected_tstep)]
	# tot_prevallist is a list of total prevalence for each time step
	tot_prevallist = [len(infected_tstep)]
	# ORlist is a list of ORs for each time step
	ORlist = []

### simulation ###
	while infected_tstep:
		tstep += 1
				
		# count to keep track of incidence
		incid_ct = 0
				
		# S to I
		suscep_tstep = [u for u in states if states[u] == 's']
		for u in suscep_tstep:
			beta = calculate_beta(T, gamma, u, dict_node_age, dict_age_susceptibility)
			
			# states[u] == 's' condition is extraneous
			if states[u] == 's' and rnd.random() < (1- np.exp(-beta*infected_neighbors(G, u, states))): 
				states[u] = 'i'
				incid_ct += 1
				# save infection time step: one col per node, one sim per row
				I_tstep_savelist[int(u)-1] = tstep

		# I to R
		for v in infected_tstep:
			# states[v] == 'i' condition is extraneous
			if states[v] == 'i' and rnd.random() < gamma:
				states[v] = 'r'
				# save recovery time step: one col per node, one sim per row
				R_tstep_savelist[int(v)-1] = tstep
		
		# update list of currently infected nodes for next tstep
		infected_tstep = [node for node in states if states[node] == 'i']
		
### metrics by time step ###
		# 1) track total incidence for each time step
		tot_incidlist.append(incid_ct)
		
		# 2) track total prevalence for each time step
		tot_prevallist.append(len(infected_tstep))
			
		# 3) return a list of ORs for each time step
		OR = calc_OR_from_list(dict_node_age, infected_tstep)
		ORlist.append(OR)

### metrics over entire simulation ###
	# 1) For list of nodes that were infected over the entire simulation, look at indexes in I_tstep_savelist or R_tstep_savelist that are not float('nan')s. Remember that index = node number - 1
	recovered = [node for node in states if states[node] == 'r']
	# 2) infected children
	rec_child_n = float(len([node for node in recovered if dict_node_age[node] == '3']))
	# 3) infected adults
	rec_adult_n = float(len([node for node in recovered if dict_node_age[node] == '4']))

	# 4) OR dict filtered to include only time points with a substantial number of infections
	min_tstep, max_tstep = filter_time_points(tot_incidlist, recovered, incl_min, incl_max)
	filtered_ORlist = [float('NaN') if (num < min_tstep or num > max_tstep) else OR for num, OR in enumerate(ORlist)]

	# 5) total OR value
	ORval_total = calc_OR_from_list(dict_node_age, recovered)


	### return data structures ###
	return rec_child_n, rec_adult_n, len(recovered), ORlist, tot_incidlist, tot_prevallist, filtered_ORlist, ORval_total, I_tstep_savelist, R_tstep_savelist

####################################################
def infected_neighbors(G, node, states):
    """ Calculate the number of infected neighbors for the node. 
	"""
    return sum([1 for node_i in G.neighbors(node) if states[node_i] == 'i']) 

####################################################
def calc_OR_from_list(dict_node_age, infected_nodelist):
	""" Calculate OR from list of infected nodes. Return OR as 1.0 if numerator or denominator are 0. 
	"""
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
def filter_time_points(tot_incidlist, recovered, incl_min, incl_max):
	""" In a time-based epidemic simulation, filter OR data points that may not be valid due to the small number of infected individuals at the beginning and end of an epidemic. For each simulation, define this period as the time period between which the percentage of infecteds is incl_min (float between 0 and 1) and incl_max (float between 0 and 1) percent of the cumulative epidemic size for that simulation. 
	"""
	
	# list of included time steps
	incl_tsteps = []
	
	for ct in xrange(len(tot_incidlist)):
		# calculate cumulative number of cases at each time step
		cum_ct = float(sum(tot_incidlist[:(ct+1)]))
		if cum_ct/len(recovered) > incl_min and cum_ct/len(recovered) < incl_max:
			incl_tsteps.append(ct)
	
	if incl_tsteps:
		return min(incl_tsteps), max(incl_tsteps)
	else:
		return 0, 0

####################################################
def define_epi_time(dict_epiincid, beta, align_proportion):
	""" For time-based epidemic simulations, identify the tstep at which the simulation reaches align_proportion proportion of cumulative infections for the epidemic. The goal is to align the epidemic trajectories at the point where the tsteps are comparable for each simulation. Returns dict where key = beta and value = list of tstep at which simulations that reached epidemic sizes attained 5% cumulative infections. 
	"""
	
	dummykeys = [k for k in dict_epiincid if k[0] == beta]
	# dict_dummyalign_tstep[beta] = [5%cum-inf_tstep_sim1, 5%cum-inf_tstep_sim2..]
	dict_dummyalign_tstep = defaultdict(list)
	
	for k in dummykeys:
		cum_infections = [sum(dict_epiincid[k][:index+1]) for index in xrange(len(dict_epiincid[k]))]
		prop_cum_infections = [float(item)/cum_infections[-1] for item in cum_infections]
		dict_dummyalign_tstep[beta].append(bisect.bisect(prop_cum_infections, align_proportion))
	
	match_tstep = int(round(np.mean(dict_dummyalign_tstep[beta])))
	
	return dict_dummyalign_tstep, match_tstep, dummykeys

####################################################
def recreate_incid(extractfile, zipname, epi_size, child_nodes, adult_nodes):
	""" For time-based epidemic simulations, recreate age-specific incidence from simulation output file of timestep at which each node got infected, where column indexes are node IDs minus one and rows are simulation results. child_nodes and adult_nodes are binary lists indicating whether the nodeID is a child/adult or not. Function returns dictionary with simnumber and age class as a tuple key and values are incidence at each tstep for a single simulation. Note: patient zero will not be counted in incidence.
	"""
	
	# calculate sizes of child and adult populations
	child_size, adult_size = float(sum(child_nodes)), float(sum(adult_nodes))
	
	# import output file from ziparchive
	with zipfile.ZipFile(zipname, 'r') as zf:
		f = zf.open(extractfile)
	
	# dict to store results data
	dict_Itstep = defaultdict(list)
	simnumber = 0
	
	# import data to dict
	for line in f:
		vallist = line.split(',')
		# dict_Itstep[simnumber] = [node1 infection tstep, node2 infection tstep...]
		# float('nan') if node was not infected
		# incorporate only sims where sim reached episize - count number of individuals infected in simulation
		if sum([0 if math.isnan(float(val)) == True else int(val) for val in vallist]) > epi_size:
			dict_Itstep[simnumber] = [0 if math.isnan(float(val)) == True else int(val) for val in vallist]
		simnumber += 1
	
	# dict to store child and adult incidence
	# dict_incid[(simnumber, 'C' or 'A')] = [C or A incid at tstep 0, C or A incid at tstep 1...]
	dict_incid = defaultdict(list)
	
	for simnum in dict_Itstep:
		# shows child and adult tsteps only - patient zero information will appear as 0
		c_tsteps = [c*t for c, t in zip(child_nodes, dict_Itstep[simnum])]
		a_tsteps = [a*t for a, t in zip(adult_nodes, dict_Itstep[simnum])]
# 		# replace NaNs with 0
# 		c_tsteps = [0 if tstep == float('nan') else tstep for tstep in c_tsteps]
# 		a_tsteps = [0 if tstep == float('nan') else tstep for tstep in a_tsteps]
		
		# create dummy lists for dict_incid values
		dummy_c_incid = [0] * np.nanmax(c_tsteps + a_tsteps)
		dummy_a_incid = [0] * np.nanmax(c_tsteps + a_tsteps)
	
		# index of dummy_c_incid is tstep and value of c_tsteps is tstep
		for t, inf in enumerate(dummy_c_incid):
			# we lost the information about age of patient zero, so just set incid for children and adults to 0 at t=0
			if t == 0:
				dummy_c_incid[t] = 0
				dummy_a_incid[t] = 0
			else:
				dummy_c_incid[t] = sum([True if int(ts) == t else False for ts in c_tsteps])/child_size * 100
				dummy_a_incid[t] = sum([True if int(ts) == t else False for ts in a_tsteps])/adult_size * 100
		
		dict_incid[(simnum, 'C')] = dummy_c_incid
		dict_incid[(simnum, 'A')] = dummy_a_incid
	
	return dict_incid


####################################################
def calculate_beta(T, gamma, susc_node, dict_node_age, dict_age_susceptibility):
	""" Calculate beta for a given susceptible node based on T and a susceptibility value for the susceptible node. Susceptibility values range from 0 (not susceptible at all) to 1 (completely susceptible). Susceptibility values are assigned by age class. This calculation operates under the assumption that T = infectivity * susceptibility.
	"""

	age_class = dict_node_age[susc_node]
	T_mod = T * dict_age_susceptibility[age_class]
	beta_mod = (-T_mod * gamma)/(T_mod - 1)
	
	# return beta value based on modified T value
	return beta_mod
	






















