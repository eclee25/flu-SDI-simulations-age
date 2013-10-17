#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/16/13
###Function: perform basic functions to explore the properties of the network

###Import data: Vancouver network files

###Command Line:  
##############################################


### notes ###
# codebook of age class codes
# '1' - Toddlers: 0-2
# '2' - Preschool: 3-4
# '3' - Children: 5-18
# '4' - Adults: 19-64
# '5' - Seniors: 65+ (community)
# '6' - Elders: 65+ (nursing home)
# There are only 94 "elders" in the Vancouver network, and they all reside in one nursing home, so they can be combined with the seniors for analysis purposes (all_elderly).

### packages/modules ###
import networkx as nx
import numpy as np

## local modules ##



##############################################
def calc_Tc(Graph):
	''' Return T_critical (R0 = 1) for a given network. Note that R0 = T * average excess degree.
	'''
	numerator = float(0)
	denominator = float(0)
	for node in Graph.nodes():
		degree = Graph.degree(node)
		numerator = numerator + (degree * (degree - 1))
		denominator = denominator + degree
	avg_excess_degree = (numerator/Graph.order())/(denominator/Graph.order())
	return 1/avg_excess_degree

##############################################
def calc_Tavg(susc, T, Graph, dict_node_age):
	''' Return average T value applied to network when susceptibility differs in adults. Recall that there will be two sets of T values for each edge; the T value for transmission to an adult node will be lower than that of all other T values.
	'''
	adult_edges = sum([len(Graph.neighbors(node)) for node in dict_node_age if dict_node_age[node] == '4'])
	network_edges = float(2) * len(Graph.edges())
	nadult_edges = network_edges - adult_edges
	T_avg = ((T * nadult_edges) + (T * susc * adult_edges)) / network_edges
	return T_avg
	











