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
	''' Time-based age-structured simulation function that returns the number of children, adults, and total infected individuals during the simulation. 
	'''
	
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	# Randomly choose one node as patient zero
	p_zero = rnd.choice(G.nodes()) 
	states[p_zero] = 'i'
	
	# define cumulative percentiles that filter OR time points
	incl_min, incl_max = 0.2, 0.8
	
	# keep track of infections over time
	# time steps begin at 0
	tstep = 0
	infected_tstep = [p_zero]
	
	# record infection timestep for patient zero
	I_tstep_savelist = [float('nan') for n in G.order()]
	I_tstep_savelist[int(p_zero)-1] = tstep # node 1 will be in column index 0 (nodes number from 1 to 10304)
	
	# create list to record recovery timesteps
	R_tstep_savelist = [float('nan') for n in G.order()]

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
	""" In a time-based epidemic simulation, filter OR data points that may not be valid due to the small number of infected individuals at the beginning and end of an epidemic. For each simulation, define this period as the time period between which the percentage of infecteds is incl_min (float between 0 and 1) and incl_max (float between 0 and 1) percent of the cumulative epidemic size for that simulation. """
	
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






























