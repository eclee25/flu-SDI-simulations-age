#Influenza epidemic in an urban population with age structure with vaccination
#Employ realistic, age-based vaccine coverage levels
#Progressively reduce vaccine efficacy
#Employ age-based vaccine efficacies/effectivenesses
#Compare AR to random antiviral treatment

# codebook of age class codes
# '1' - Toddlers: 0-2
# '2' - Preschool: 3-4
# '3' - Children: 5-18
# '4' - Adults: 19-64
# '5' - Seniors: 65+ (community)
# '6' - Elders: 65+ (nursing home)
# There are only 94 "elders" in the Vancouver network, and they all reside in one nursing home, so they can be combined with the seniors for analysis purposes (all_elderly).

###### ECL changes ######
# 7/22/13 ECL: descriptive comments
# 7/24/13 ECL: calculateOR function
# 8/5/13 ECL: start of major code reorganization
###########################################


### import packages and modules ###
from networkx import *
from random import *
from pylab import mean, hist, show
from collections import defaultdict
import random as rnd
import sys
import pretty_print as gen



### local functions ###
# 1) percolation function
def percolate(G,i,s):
    # start of simulation
    states = dict([(node, 's') for node in G.nodes()]) # all nodes are assigned 's' at start
    p_zero = choice(G.nodes()) # pick a random node as patient zero
    states[p_zero] = 'i' # assign infected state to patient zero
    
    # keep track of nodes that are infected and recovered during sim
    infected, recovered = [p_zero],[]
    
    # keep track of nodes that have been infected by age group during sim
    toddlers, preschool, children = [],[],[]
    adults, seniors, elders, all_elderly = [],[],[],[]

    # keep track of vax nodes that have been infected by age group during sim
    protected = []
    p_toddlers, p_preschool, p_children = [],[],[]
    p_adults, p_seniors, p_elders, p_all_elderly = [],[],[],[]

    while len(infected) > 0:
        v = infected.pop(0)
        i = infect[v] # transmissibility value for the scenario, not age-specific
        for u in G.neighbors(v):
            s = suscept[u]  # susceptiblity value for the age group
            if states[u] == 's' and random() < (i*s):
                states[u] = 'i'
                infected.append(u)
        states[v] = 'r'
        recovered.append(v)
        if vaccine[v] == 1:
            protected.append(v)
    for node in recovered:
        if ages[node] == '1':
            toddlers.append(node)
        elif ages[node] == '2':
            preschool.append(node)
        elif ages[node] == '3':
            children.append(node)
        elif ages[node] == '4':
            adults.append(node)
        elif ages[node] == '5':
            seniors.append(node)
            all_elderly.append(node)
        else:
            elders.append(node)
            all_elderly.append(node)
    for node in protected: # Vaccinated individuals are still somewhat susceptible ECL 8/5/13
        if ages[node] == '1':
            p_toddlers.append(node)
        elif ages[node] == '2':
            p_preschool.append(node)
        elif ages[node] == '3':
            p_children.append(node)
        elif ages[node] == '4':
            p_adults.append(node)
        elif ages[node] == '5':
            p_seniors.append(node)
            p_all_elderly.append(node)
        else:
            p_elders.append(node)
            p_all_elderly.append(node)
    return (len(recovered), len(toddlers), len(preschool), len(children),
            len(adults), len(seniors), len(elders), len(all_elderly),
            len(protected), len(p_toddlers), len(p_preschool), len(p_children),
            len(p_adults), len(p_seniors), len(p_elders), len(p_all_elderly))

# Define OR calculation function
def calculateOR(child_AR, adult_AR): # defined by ECL 7/23/13
    oddsratio = (child_AR/(1-child_AR))/(adult_AR/(1-adult_AR))
    return oddsratio


### data structures ###
d_node_age = {} # key = node number, value = age class code
d_age_nodelist = defaultdict(list) # key = age class, value = list of nodes belonging to age class

### data import ###
## import Vancouver network
file = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_edges_Sarah.csv')
G = Graph()
for edge in file:
    G.add_edge(*edge.strip().split(','))
net_size = G.order()

##  import age structure
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/Data/urban_ages_Sarah.csv')
# file has node number and age category

for line in file2:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        d_node_age[node] = age

for node in d_node_age:
    if d_node_age[node] == '1':
        d_age_nodelist['toddlers'].append(node)
    elif d_node_age[node] == '2':
        d_age_nodelist['preschoolers'].append(node)
    elif d_node_age[node] == '3':
        d_age_nodelist['children'].append(node)
    elif d_node_age[node] == '4':
        d_age_nodelist['adults'].append(node)
    elif d_node_age[node] == '5':
        d_age_nodelist['seniors'].append(node) 
        d_age_nodelist['all_elderly'].append(node)  
    elif d_node_age[node] == '6':
        d_age_nodelist['elders'].append(node)
        d_age_nodelist['all_elderly'].append(node)





#Run simulation for each T-value
# tau = transmissibility
# eff = vaccine efficacy

#coverage = {0.0643: 0.20, 0.0767: 0.40}
#efficacy = {0.0643: 0.65, 0.0767: 0.65}
#prob_uninfected = {0.0643: 0.7815823810, 0.0767: 0.5476333532}

simulated_VE = {}

reduction_list = [0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00]         
# as the value in reduction_list increases, susceptibility increases
#reduction_list = [0.00,0.25,0.50,0.75,1.00]
#reduction_list = [0.05,0.95]
#reduction_list = [1.0,0.00]

numsims = 500 # number of simulations
age_list = ['toddlers','preschool','children','adults','elderly']
for reduction in reduction_list:
    #scenarios = ['seasonal']
    scenarios = ['seasonal','high eff/low cov','low eff/high cov']
    while len(scenarios) > 0:
        q = scenarios.pop(0)
        # define the parameters for each scenario
        if q == 'seasonal': #Seasonal
            tau = 0.0643
            #tau = 0.0572
#             u = 0.7815823810 # used to find the susceptibility values from the sigma calculation
            eff = 0.645
            cov = 'Seasonal'
            #Age-based coverage levels
            #overall coverage = 0.245
            cov_toddlers = 0.30
            cov_preschool = 0.20
            cov_children = 0.20
            cov_adults = 0.20
            cov_seniors = 0.65
            cov_elders = 0.90
            #Age-based susceptibility values
            #Pandemic values from sigma calculation
            sigma_dict_baseline = {'toddlers':0.52,'preschool':0.42,'children':0.33,
                                   'adults':0.28,'elderly':0.55}
            sigma_dict = {}
            for age in age_list:
                sigma_dict[age] = sigma_dict_baseline[age] + reduction*(1-sigma_dict_baseline[age])
            #print reduction, sigma_dict
        elif q == 'high eff/low cov': #Pandemic - low coverage, high efficacy
            tau = 0.0767
            #tau = 0.1132
#             u = 0.5476333532 # used to find the susceptibility values from the sigma calculation
            eff = 0.795
            cov = 'Low'
            #Age-based coverage levels
            #overall coverage = 0.290
            cov_toddlers = 0.40
            cov_preschool = 0.40
            cov_children = 0.40
            cov_adults = 0.22
            cov_seniors = 0.30
            cov_elders = 0.90
            #Age-based susceptibility values
            sigma_dict_baseline = {'toddlers':0.36,'preschool':0.27,'children':0.21,
                                   'adults':0.13,'elderly':0.36}
            sigma_dict = {}
            for age in age_list:
                sigma_dict[age] = sigma_dict_baseline[age] + reduction*(1-sigma_dict_baseline[age])
            #print reduction, sigma_dict
        elif q == 'low eff/high cov': #Pandemic - high coverage, low efficacy
            tau = 0.0767
            #tau = 0.1132
#             u = 0.5476333532 # used to find the susceptibility values from the sigma calculation
            eff = 0.645
            cov = 'High'
            #Age-based coverage levels
            #overall coverage = 0.455
            cov_toddlers = 0.60
            cov_preschool = 0.40
            cov_children = 0.40
            cov_adults = 0.40
            cov_seniors = 0.90
            cov_elders = 0.95
            #Age-based susceptibility values
            sigma_dict_baseline = {'toddlers':0.52,'preschool':0.42,'children':0.33,
                                   'adults':0.28,'elderly':0.55}
            sigma_dict = {}
            for age in age_list:
                sigma_dict[age] = sigma_dict_baseline[age] + reduction*(1-sigma_dict_baseline[age])
            #print reduction, sigma_dict

        vaccine = {} # dict indicating whether individual is vaccinated or not in the simulation
        
        #REALISTIC
        #<agegroup>_count is the list of nodes that corresponds to that age group
        number_toddlers_to_vaccinate = int(round(len(toddlers_count)*cov_toddlers))
        number_preschool_to_vaccinate = int(round(len(preschool_count)*cov_preschool))
        number_children_to_vaccinate = int(round(len(children_count)*cov_children))
        number_adults_to_vaccinate = int(round(len(adults_count)*cov_adults))
        number_seniors_to_vaccinate = int(round(len(seniors_count)*cov_seniors))
        number_elders_to_vaccinate = int(round(len(elders_count)*cov_elders))

        '''
        print number_toddlers_to_vaccinate
        print number_preschool_to_vaccinate
        print number_children_to_vaccinate
        print number_adults_to_vaccinate
        print number_seniors_to_vaccinate
        print number_elders_to_vaccinate
        '''

        todvax_list = rnd.sample(toddlers_count,number_toddlers_to_vaccinate)
        prevax_list = rnd.sample(preschool_count,number_preschool_to_vaccinate)
        chivax_list = rnd.sample(children_count,number_children_to_vaccinate)
        aduvax_list = rnd.sample(adults_count,number_adults_to_vaccinate)
        senvax_list = rnd.sample(seniors_count,number_seniors_to_vaccinate)
        eldvax_list = rnd.sample(elders_count,number_elders_to_vaccinate)

        vacc_list = (todvax_list + prevax_list + chivax_list + aduvax_list +
                     senvax_list + eldvax_list) # list of vaccinated nodes
        set_vacc = set(vacc_list)
        set_all = set(G.nodes())
        set_unvacc = set_all.difference(set_vacc)
        unvacc_list = list(set_unvacc)
        
        for node in vacc_list: # shows key and 1/0 for whether node is vaccinated or unvacc
            vaccine[node] = 1
        for node in unvacc_list:
            vaccine[node] = 0
        
        print 'Fraction Vaccinated in Network:',float(len(vacc_list))/net_size
        
        infect = {}
        for node in G.nodes():
            infect[node] = float(tau)

        suscept = {}
        # sigma_dict_baseline is different for each scenario
        # susceptibility is lower than 1 for vacc individuals only
        # as the value in reduction_list increases, susceptibility increases
        for node in G.nodes():
            if vaccine[node] == 1: # only nodes that are vaccinated have a non-one susceptibility
                if ages[node] == '1':
                    suscept[node] = float(sigma_dict['toddlers'])
                elif ages[node] == '2':
                    suscept[node] = float(sigma_dict['preschool'])
                elif ages[node] == '3':
                    suscept[node] = float(sigma_dict['children'])
                elif ages[node] == '4':
                    suscept[node] = float(sigma_dict['adults'])
                elif ages[node] == '5' or ages[node] == '6':
                    suscept[node] = float(sigma_dict['elderly'])
            else:
                suscept[node] = float(1) # susceptibility is 1 for nodes that are unvacc
        
        #Set infectivity and susceptibility
        i = 0
        s = 0

        #Run simulation
        results = []
        r_toddlers = []
        r_preschool = []
        r_children = []
        r_adults = []
        r_seniors = []
        r_elders = []
        r_all_elderly = []

        epidemics = []
        e_toddlers = []
        e_preschool = []
        e_children = []
        e_adults = []
        e_seniors = []
        e_elders = []
        e_all_elderly = []

        ar_vax = []
        #ar_age = []
        ar_vax_toddlers = []
        ar_vax_preschool = []
        ar_vax_children = []
        ar_vax_adults = []
        ar_vax_allelderly = []

        for i in range(numsims):
            z,a,b,c,d,e,f,g,y,h,j,k,l,m,n,o = percolate(G,i,s)
            results.append(z)
            r_toddlers.append(a)
            r_preschool.append(b)
            r_children.append(c)
            r_adults.append(d)
            r_seniors.append(e)
            r_elders.append(f)
            r_all_elderly.append(g)
            if z > 515:
                epidemics.append(z)
                e_toddlers.append(a)
                e_preschool.append(b)
                e_children.append(c)
                e_adults.append(d)
                e_seniors.append(e)
                e_elders.append(f)
                e_all_elderly.append(g)
                
                ar_vax.append(y)
                ar_vax_toddlers.append(h)
                ar_vax_preschool.append(j)
                ar_vax_children.append(k)
                ar_vax_adults.append(l)
                ar_vax_allelderly.append(o)
                #ar_age.append(m)

        #View results without saving to text files
        print 'Resulted in epidemic (%):', len(epidemics)/float(numsims)*100 # percentage of simulations that resulted in epidemics
        print 'T:',tau
        print 'Coverage:',cov
        print 'Reduction in vax efficacy:',reduction
        print
        
        if len(epidemics) > 0:
            print 'Mean epidemic size:', float(mean(epidemics))/net_size
            '''
            if len(e_toddlers) > 0:
                print 'Toddlers:', float(mean(e_toddlers))/len(toddlers_count)
            else:
                print 'Toddlers: 0'
            if len(e_preschool) > 0:
                print 'Preschool:', float(mean(e_preschool))/len(preschool_count)
            else:
                print 'Preschool: 0'
            if len(e_children) > 0:
                print 'Children:', float(mean(e_children))/len(children_count)
            else:
                print 'Children: 0'
            if len(e_adults) > 0:
                print 'Adults:', float(mean(e_adults))/len(adults_count)
            else:
                print 'Adults: 0'
            if len(e_seniors) > 0:
                print 'Seniors:', float(mean(e_seniors))/len(seniors_count)
            else:
                print 'Seniors: 0'
            if len(e_elders):
                print 'Elders:', float(mean(e_elders))/len(elders_count)
            else:
                print 'Elders: 0'
            '''  
        else:
            print 'No epidemics'
        #print 'Probability of epidemic:', float(len(epidemics))/float(len(results))
        #print 'Size x Prob:', (float(mean(epidemics))/net_size)*(float(len(epidemics))/float(len(results)))
        
        
        '''
        #CHECK
        #print 'CHECK'
        #print 'Vaccinated:',float(len(vacc_list))/net_size
        #print
        #print 'Expected: 0.071/0.082/0.142'
        #print 'AR total:',float(mean(epidemics))/net_size
        '''

        #EFFICACY ANALYSIS
        #print 'Total:'
        ar_vacc_total = float(mean(ar_vax))/len(vacc_list)
        ar_unvacc_total = (float(mean(epidemics))-float(mean(ar_vax)))/(net_size-len(vacc_list))
        #print ar_vacc_total
        #print ar_unvacc_total
        calc_VE_total = ((ar_unvacc_total - ar_vacc_total)/ar_unvacc_total)*100
        #print calc_VE_total
        #print
        #print 'By age:'
        #print '    Toddlers:'
        ar_vacc_toddlers = float(mean(ar_vax_toddlers))/len(todvax_list)
        ar_unvacc_toddlers = (float(mean(e_toddlers))-float(mean(ar_vax_toddlers)))/(len(toddlers_count)-len(todvax_list))
        #print '    ',ar_vacc_toddlers
        #print '    ',ar_unvacc_toddlers
        calc_VE_toddlers = ((ar_unvacc_toddlers - ar_vacc_toddlers)/ar_unvacc_toddlers)*100
        #print '    ',calc_VE_toddlers
        #print
        #Calculate for all other age groups, then print at the end
        ar_vacc_preschool = float(mean(ar_vax_preschool))/len(prevax_list)
        ar_unvacc_preschool = (float(mean(e_preschool))-float(mean(ar_vax_preschool)))/(len(preschool_count)-len(prevax_list))
        ar_vacc_children = float(mean(ar_vax_children))/len(chivax_list)
        ar_unvacc_children = (float(mean(e_children))-float(mean(ar_vax_children)))/(len(children_count)-len(chivax_list))
        ar_vacc_adults = float(mean(ar_vax_adults))/len(aduvax_list)
        ar_unvacc_adults = (float(mean(e_adults))-float(mean(ar_vax_adults)))/(len(adults_count)-len(aduvax_list))
        ar_vacc_allelderly = float(mean(ar_vax_allelderly))/(len(senvax_list)+len(eldvax_list))
        ar_unvacc_allelderly = (float(mean(e_all_elderly))-float(mean(ar_vax_allelderly)))/(len(all_elderly_count)-(len(senvax_list)+len(eldvax_list)))
        
        calc_VE_preschool = ((ar_unvacc_preschool - ar_vacc_preschool)/ar_unvacc_preschool)*100
        calc_VE_children = ((ar_unvacc_children - ar_vacc_children)/ar_unvacc_children)*100
        calc_VE_adults = ((ar_unvacc_adults - ar_vacc_adults)/ar_unvacc_adults)*100
        calc_VE_allelderly = ((ar_unvacc_allelderly - ar_vacc_allelderly)/ar_unvacc_allelderly)*100

        '''
        print '    Preschool:'
        #print '    ',ar_vacc_preschool
        #print '    ',ar_unvacc_preschool
        print '    ',calc_VE_preschool
        #print
        print '    Children:'
        #print '    ',ar_vacc_children
        #print '    ',ar_unvacc_children
        print '    ',calc_VE_children
        #print
        print '    Adults:'
        #print '    ',ar_vacc_adults
        #print '    ',ar_unvacc_adults
        print '    ',calc_VE_adults
        #print
        print '    All Elderly:'
        #print '    ',ar_vacc_allelderly
        #print '    ',ar_unvacc_allelderly
        print '    ',calc_VE_allelderly
        #print
        '''
        
        #print
        #print 'Average:'
        mean_calc_VE = (len(toddlers_count)*calc_VE_toddlers+len(preschool_count)*calc_VE_preschool+len(children_count)*calc_VE_children+
                        len(adults_count)*calc_VE_adults+len(all_elderly_count)*calc_VE_allelderly)/net_size
        print "overall vax effec. (%):", mean_calc_VE
        simulated_VE[(q,reduction)] = mean_calc_VE
        print '--------------------'
        print
        

        #Print results to text files
        '''
        l = results
        filename = 'realistic_TIV_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        '''
        '''
        l = r_toddlers
        filename = 'realistic_TIV_toddlers_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = r_preschool
        filename = 'realistic_TIV_preschool_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = r_children
        filename = 'realistic_TIV_children_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = r_adults
        filename = 'realistic_TIV_adults_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = r_seniors
        filename = 'realistic_TIV_seniors_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = r_elders
        filename = 'realistic_TIV_elders_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = r_all_elderly
        filename = 'realistic_TIV_allelderly_T%r_%s_%s.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        
        l = low_symptom
        filename = 'realistic_TIV_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        
        l = ls_toddlers
        filename = 'realistic_TIV_toddlers_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = ls_preschool
        filename = 'realistic_TIV_preschool_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = ls_children
        filename = 'realistic_TIV_children_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = ls_adults
        filename = 'realistic_TIV_adults_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = ls_seniors
        filename = 'realistic_TIV_seniors_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = ls_elders
        filename = 'realistic_TIV_elders_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        l = ls_all_elderly
        filename = 'realistic_TIV_allelderly_T%r_%s_%s_ls.txt'%(tau,cov,reduction)
        gen.print_list_to_file(l,filename)
        '''

'''
scenarios_2 = ['seasonal','high eff/low cov','low eff/high cov']
for scenario in scenarios_2:
    print scenario
    for reduction in reduction_list:
        print reduction, ':', simulated_VE[(scenario,reduction)]
    print
'''

#Epidemic = more than 515 people infected
