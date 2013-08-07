#Influenza epidemic in an urban population with age structure with vaccination
#Employ realistic, age-based vaccine coverage levels
#Progressively reduce vaccine efficacy
#Employ age-based vaccine efficacies/effectivenesses
#Compare AR to random antiviral treatment

# Toddlers: 0-2
# Preschool: 3-4
# Children: 5-18
# Adults: 19-64
# Seniors: 65+ (community)
# Elders: 65+ (nursing home)
# There are only 94 "elders" in the Vancouver network, and they all reside in one nursing home, so they can generally be combined with the seniors for analysis purposes (all_elderly).



from networkx import *
from random import *
from pylab import mean, hist, show
import random as rnd
import sys
# sys.path.append('../../../../Models') # commented out ECL 7/22
import pretty_print as gen

#Add in age structure
file2 = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/urban_ages_Sarah.csv')
ages = {}
for line in file2:
    new_line = line.split()
    for line in new_line:
        node, age = line.split(',')
        ages[node] = age

#Form lists of nodes in each age group
toddlers_count = []
preschool_count = []
children_count = []
adults_count = []
seniors_count = []
elders_count = []
all_elderly_count = []

for node in ages:
    if ages[node] == '1':
        toddlers_count.append(node)
    elif ages[node] == '2':
        preschool_count.append(node)
    elif ages[node] == '3':
        children_count.append(node)
    elif ages[node] == '4':
        adults_count.append(node)
    elif ages[node] == '5':
        seniors_count.append(node)
    else:
        elders_count.append(node)
all_elderly_count = seniors_count + elders_count

#Insert network
file = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/Age_Based_Simulations/urban_edges_Sarah.csv')
G = Graph()
for edge in file:
    G.add_edge(*edge.strip().split(','))
net_size = G.order()

#Define percolation function
def percolate(G,i,s):
    states = dict([(node, 's') for node in G.nodes()])

    p_zero = choice(G.nodes())
    states[p_zero] = 'i'
    infected = [p_zero]
    recovered = []

    toddlers = []
    preschool = []
    children = []
    adults = []
    seniors = []
    elders = []
    all_elderly = []

    protected = []
    p_toddlers = []
    p_preschool = []
    p_children = []
    p_adults = []
    p_seniors = []
    p_elders = []
    p_all_elderly = []

    while len(infected) > 0:
        v = infected.pop(0)
        i = infect[v] # transmissibility value for the scenario, not age-specific
        for u in G.neighbors(v):
            s = suscept[u] # susceptiblity value for the age group
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
    for node in protected:
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

#Run simulation for each T-value
# tau = transmissibility
# eff = vaccine efficacy

age_list = ['toddlers','preschool','children','adults','elderly']
scenarios = ['seasonal','high eff/low cov','low eff/high cov']
while len(scenarios) > 0:
    q = scenarios.pop(0)
    if q == 'seasonal': #Seasonal
        tau = 0.0643
        #tau = 0.0572
        eff = 0.645 # reference value of weighted average of vax efficacy
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
        sigma_dict = {'toddlers':0.52,'preschool':0.42,'children':0.33,
                      'adults':0.28,'elderly':0.55}
    elif q == 'high eff/low cov': #Pandemic - low coverage, high efficacy
        tau = 0.0767
        #tau = 0.1132
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
        sigma_dict = {'toddlers':0.36,'preschool':0.27,'children':0.21,
                      'adults':0.13,'elderly':0.36}
    elif q == 'low eff/high cov': #Pandemic - high coverage, low efficacy
        tau = 0.0767
        #tau = 0.1132
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
        sigma_dict = {'toddlers':0.52,'preschool':0.42,'children':0.33,
                      'adults':0.28,'elderly':0.55}
        
    vaccine = {}
    
    #REALISTIC
    #import random as rnd
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
                 senvax_list + eldvax_list)
    set_vacc = set(vacc_list)
    set_all = set(G.nodes())
    set_unvacc = set_all.difference(set_vacc)
    unvacc_list = list(set_unvacc)
    
    for node in vacc_list:
        vaccine[node] = 1 # vaccine dict designates the vax and non-vax nodes
    for node in unvacc_list:
        vaccine[node] = 0
    
    #print 'Vaccinated:',float(len(vacc_list))/net_size
    
    infect = {}
    for node in G.nodes():
        infect[node] = float(tau)

    suscept = {}
    for node in G.nodes():
        if vaccine[node] == 1:
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
            suscept[node] = float(1)
    
    #Set infectivity and susceptibility
    i = 0 # important to clear the values from the previous iterations
    s = 0

    #Run simulation
    results, r_toddlers, r_preschool, r_children, r_adults, r_seniors, r_elders, r_all_elderly = [],[],[],[],[],[],[],[] # number of recovered for all simulations

    epidemics, e_toddlers, e_preschool, e_children, e_adults, e_seniors, e_elders, e_all_elderly = [],[],[],[],[],[],[],[] # number of recovered for simulations that become epidemics (number of infecteds is over 515)

    
	for i in range(500):
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


    #View results without saving to text files
    print 'T:',tau
    print 'Coverage:',cov
  
    if len(epidemics) > 0:
        print 'Mean epidemic size:', float(mean(epidemics))/net_size
        
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
        child_AR_e = [float(val)/len(children_count) for val in e_children] # child attack rate for each simulation that reaches epidemic proportions
        adult_AR_e = [float(val)/len(adults_count) for val in e_adults] # adult attack rate for each simulation that reaches epidemic proportions
        print "number of children infected during epidemic", e_children, len(children_count)
        print "number of adults infected during epidemic", e_adults, len(adults_count)
        #print "child attackrates", child_AR
        #print "adult attackrates", adult_AR
        OR = map(calculateOR, child_AR_e, adult_AR_e)
        print "ORs (epidemics only):", OR	
        print 'mean OR (epidemics only):', mean(OR) 
    else:
        print 'No epidemics'
    #print 'Probability of epidemic:', float(len(epidemics))/float(len(results))
    #print 'Size x Prob:', (float(mean(epidemics))/net_size)*(float(len(epidemics))/float(len(results)))
    print
    
'''
    #CHECK
    print
	print 'CHECK'
    print 'Vaccinated:',float(len(vacc_list))/net_size
    print
    print 'Expected: 0.071/0.082/0.142'
    print 'AR total:',float(mean(epidemics))/net_size
'''  

    #Print results to text files
'''
    l = results
    filename = 'realistic_TIV_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)


    l = r_toddlers
    filename = 'realistic_TIV_toddlers_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)
    l = r_preschool
    filename = 'realistic_TIV_preschool_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)
    l = r_children
    filename = 'realistic_TIV_children_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)
    l = r_adults
    filename = 'realistic_TIV_adults_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)
    l = r_seniors
    filename = 'realistic_TIV_seniors_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)
    l = r_elders
    filename = 'realistic_TIV_elders_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)
    l = r_all_elderly
    filename = 'realistic_TIV_allelderly_T%r_%s.txt'%(tau,cov)
    gen.print_list_to_file(l,filename)
'''    


#Epidemic = more than 515 people infected
