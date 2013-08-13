#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 8/7/13
###Purpose:
#### number of child and adult nodes in network
#### separate epidemics from results
#### calculate incidence of children and adults??
#### random vax strategy
#### vax efficacy perc function for age-structured networks
#### gen perc function for age-structured networks

###Import data: 

###Command Line: python 
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


# functions
####################################################
# number of child and adult nodes in network
def child_adult_size(dict_node_age):
	child_size = len([node for node in dict_node_age if dict_node_age[node] =='3'])
	adult_size = len([node for node in dict_node_age if dict_node_age[node] =='4'])
	return child_size, adult_size

###################################################
# separate epidemics from all results
def epidemicsonly(dict_simresults, episize):
	return dict((k, dict_simresults[k]) for k in dict_simresults.keys() if dict_simresults[k][2] > episize)

# ###################################################
# # calculate incidence of children and adults 
# # won't return anything, dictionary items are added to an empty dictionary -- is this okay for readability? (not currently used in the main code)
# def child_adult_incid(dict_simepi, dict_epiincid, child_size, adult_size):
# 	for key in dict_simepi:
# 		dict_epiincid[(key[0], 'C')].append(dict_simepi[key][0]/child_size)
# 		dict_epiincid[(key[0], 'A')].append(dict_simepi[key][1]/adult_size)

####################################################
# random vax strategy with vaxcov% coverage and 100% efficacy
def rnd_vax(G, vaxcov):
	for n in G.nodes():
		if rnd.random() < vaxcov:
			G.remove_edges_from(G.edges(n))

####################################################
# age structured percolation function with random vax strategy with vaxcov% coverage and defined susceptibility multiplier, where 1 is fully susceptible
def perc_age_vaxeff(G, dict_node_age, T, vaxcov_scen, susceptibility):

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
	while len(infected) > 0:
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
	
	# infected children
	rec_child_n = float(len([node for node in recovered if dict_node_age[node] == '3']))
	# infected adults
	rec_adult_n = float(len([node for node in recovered if dict_node_age[node] == '4']) )
	
	# number recovered children, adults, and total pop, and number of initially vaccinated individuals
	return rec_child_n, rec_adult_n, len(recovered), len(vaccinated)

####################################################
# general age-structured percolation function 
def perc_age_gen(G, dict_node_age, T):
	# set initial conditions
	states = dict([(node, 's') for node in G.nodes()])
	p_zero = rnd.choice(G.nodes()) # Randomly choose one node as patient zero
	states[p_zero] = 'i'
	infected    = [p_zero] 
	recovered   = [] 
	
	# simulation
	while len(infected) > 0:
		v = infected.pop(0)
		for u in G.neighbors(v):
			if states[u] == 's' and rnd.random() <T: # make sure node u is still susceptible
				states[u] = 'i'
				infected.append(u)
		states[v] = 'r'
		recovered.append(v)
	
	# infected children
	rec_child_n = float(len([node for node in recovered if dict_node_age[node] == '3']))
	# infected adults
	rec_adult_n = float(len([node for node in recovered if dict_node_age[node] == '4']) )
	
	# number recovered in children and adults and total
	return rec_child_n, rec_adult_n, len(recovered)







### Inefficient percolations - do NOT use ###
# #############################
# # DICT ONLY: general age-structured percolation function
# def perc_age_gen(G, dict_node_age, T):
# 	# set initial conditions
# 	states = dict([(node, 's') for node in G.nodes()])
# 	p_zero = rnd.choice(G.nodes()) # Randomly choose one node as patient zero
# 	states[p_zero] = 'i'
# 	
# 	# simulation
# 	while len([states[node] for node in states if states[node]=='i']) > 0:
# 		v = rnd.choice([node for node in states if states[node] == 'i'])
# 		for u in G.neighbors(v):
# 			if states[u] == 's' and rnd.random() <T: # make sure node u is still susceptible
# 				states[u] = 'i'
# 		states[v] = 'r'
# 	
# 	# infected children
# 	rec_child_n = float(len([node for node in states if dict_node_age[node] == '3' and states[node] == 'r']))
# 	# infected adults
# 	rec_adult_n = float(len([node for node in states if dict_node_age[node] == '4' and states[node] == 'r']))
# 	
# 	# return incidence in children and adults
# 	return rec_child_n, rec_adult_n, len([node for node in states if states[node]=='r'])

##################################
## LIST ONLY: general age-structured percolation function
## this method was magnitudes slower than the dictionary method of percolation
# def perc_age_gen(G, dict_node_age, T):
# 	# set initial conditions
# 	p_zero = rnd.choice(G.nodes()) # Randomly choose one node as patient zero
# 	infected    = [p_zero]
# 	recovered   = []
# 	
# 	# simulation
# 	while len(infected) > 0:
# 		v = infected.pop(0)
# 		for u in G.neighbors(v):
# 			# We need to make sure Node u is still susceptible
# 			if u not in infected and u not in recovered and rnd.random()<T:
# 				infected.append(u)
# 		recovered.append(v)
# 	
# 	# infected children
# 	rec_child_n = float(len([node for node in recovered if dict_node_age[node] == '3']))
# 	# infected adults
# 	rec_adult_n = float(len([node for node in recovered if dict_node_age[node] == '4']) )
# 	
# 	# return incidence in children and adults
# 	return rec_child_n, rec_adult_n, len(recovered)
