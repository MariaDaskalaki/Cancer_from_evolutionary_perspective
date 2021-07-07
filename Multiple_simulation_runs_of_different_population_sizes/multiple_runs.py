#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:40:32 2020

@author: maria
"""
################# Code for plotting fixation time boxplots and distributions and invasion time boxplots and distribution for different population sizes #########################

import csv
from itertools import chain
from itertools import accumulate
import re
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import numpy as np

with open ("/home/daskalaki/synology/daskalaki/population_vs_fixation_time.txt", "w+") as final_file:
        final_file.write("Population_size\tMean_fixation_time\tProbability_of_fixation\tMean_invasion_fixation_time\tInvasion_probability\n")

nruns=100
path_list=["/home/daskalaki/synology/daskalaki/100pop_new", "/home/daskalaki/synology/daskalaki/300pop_new", "/home/daskalaki/synology/daskalaki/500pop_new", "/home/daskalaki/synology/daskalaki/800pop_new", "/home/daskalaki/synology/daskalaki/1000pop_new"]
default_dictionary=defaultdict(list)
default_cancerous_dictionary=defaultdict(list)

for path in path_list:
    fixation_list=[]
    cancerous_fixation_list=[]
    path_for_population_size= path+"/results1/Number of mutants.txt"
    
    with open (path_for_population_size, "r") as f_pop:
        data=f_pop.readlines()
        initial_population_size=int(re.sub("[^0-9]", "", data[0])) # initial population size
        print("The population size is {}".format(initial_population_size))


    for i in range(nruns):
        path_for_run= path+"/results{}/Frequencies_CA.txt".format(i+1)
        path_for_run_cancerous=path+"/results{}/Number of mutants.txt".format(i+1)
        generation_count=0
        with open (path_for_run, "r") as f:
            data=f.readlines()
            generation_list=[int(line.split()[1]) for line in data]
            fixation_generation=generation_list[-1] # in which generation fixation has hapenned
            #fixation_list.append(fixation_generation)
            genot_freq_list=[line.split()[2:] for line in data]
            #print(genot_freq_list)
            #print("There has been fixation in generation ",count)
            genot_freq_list_final=list(chain.from_iterable(genot_freq_list)) # list of frequencies for each genotype
            #print(genot_freq_list_final)
            genot_list=[i.split(':')[0] for i in genot_freq_list_final]
            #print (genot_list)
            freq_list=[float(i.split(':')[1]) for i in genot_freq_list_final]
            #print(freq_list)
            split_length=[len(i) for i in genot_freq_list]
            #print(split_length)
            freq_list_new=[freq_list[x - y: x] for x, y in zip(accumulate(split_length), split_length)]
            #print(freq_list_new)
            fixation_list.append(fixation_generation) # list of fixation times for each simulation run of each population size
            key=initial_population_size
            value=fixation_generation
            default_dictionary[key].append(value)
            
    
        with open (path_for_run_cancerous, "r") as file:
                main_data=file.readlines()[1:]
                for line in main_data:
                    find=False
                    if int(line.split()[2]) == initial_population_size: # stands for invasion of cancerous cells in population
                        find=True
                    if find==True:
                        fixed_generation=int(line.split()[1])
                        cancerous_fixation_list.append(fixed_generation)
                        key_cancerous=initial_population_size
                        value_cancerous=fixed_generation
                        default_cancerous_dictionary[key_cancerous].append(value_cancerous)
                        print("The fixation of cancerous has happened in generation {}".format(fixed_generation))
                        break
                if find==False:
                    print ("The has been no fixation of cancerous")
    print (cancerous_fixation_list) # list of invasion time for each simulation run of each population size
    print ("There have been {} fixation of cancerous".format(len(cancerous_fixation_list)))
    mean_cancerous_fixation_time=sum(cancerous_fixation_list)/len(cancerous_fixation_list)
    print ("The mean fixation time of cancerous for population size {} is {}".format(initial_population_size,mean_cancerous_fixation_time))
    invasion_probability=len(cancerous_fixation_list)/nruns # mean invasion time for each population size
    print ("The invasion probability for population size {} is {}".format(initial_population_size, invasion_probability))


    print(fixation_list)
    #print("There are {} fixations".format(fixation_count))
    print("There are {} fixations".format(len(fixation_list)))

    mean_fixation_time=sum(fixation_list)/len(fixation_list)
    probability_of_fixation=len(fixation_list)/nruns # mean fixation time for each population size
    print("The mean fixation time is {}".format(mean_fixation_time))
    print("The probability of fixation is {}".format(probability_of_fixation))
    
   

    with open ("/home/daskalaki/synology/daskalaki/population_vs_fixation_time.txt", "a+") as final_file:
        #final_file.write("Population_size\tMean_fixation_time\n")
        final_file.write("{}\t{}\t{}\t{}\t{}\n".format(initial_population_size, mean_fixation_time, probability_of_fixation, mean_cancerous_fixation_time, invasion_probability))
        
print (default_dictionary)
mean_values=[np.mean(i) for i in default_dictionary.values()]
print(mean_values)

print ("The default_cancerous_dictionary is")
print (default_cancerous_dictionary)
mean_cancerous_values=[np.mean(i) for i in default_cancerous_dictionary.values()]
print ("The cancerous mean values are")
print (mean_cancerous_values)


fig, ax = plt.subplots()
box=ax.boxplot(default_dictionary.values())
ax.set_xticklabels(default_dictionary.keys())
plt.plot([1,2,3,4,5], mean_values, c='b', lw=2)
plt.xlabel("Population sizes")
plt.ylabel("Fixation times")
plt.title("Fixation time boxplots- Mean fixation time distribution")
plt.savefig("im1_new.pdf")
plt.show()

fig, ax = plt.subplots()
box=ax.boxplot(default_cancerous_dictionary.values())
ax.set_xticklabels(default_cancerous_dictionary.keys())
plt.plot([1,2,3,4,5], mean_cancerous_values, c='b', lw=2)
plt.xlabel("Population sizes")
plt.ylabel("Invasion times")
plt.title("Invasion time boxplots- Mean invasion fixation of cancerous time distribution")
plt.savefig("im2_new.pdf")
plt.show()



data_final=pd.read_csv("/home/daskalaki/synology/daskalaki/population_vs_fixation_time.txt", sep='\t', lineterminator='\n', header=0)
#print(data_final)
data_final.plot(x='Population_size', y='Mean_fixation_time')
plt.ylabel("Generations")
plt.title("beta_gom_time_dep_stochastic_nonherited model")
plt.savefig("im3_new.pdf")
plt.show()

data_final.plot(x='Population_size', y='Probability_of_fixation')
plt.ylabel("Fixation Probabilities")
plt.title("Fixation probabilities distribution in beta_gom_time_dep_stochastic_nonherited_model")
plt.savefig("im4_new.pdf")
plt.show()

data_final.plot(x='Population_size', y='Mean_invasion_fixation_time')
plt.ylabel("Generations")
plt.title("Mean invasion time in beta_gom_time_dep_stochastic_nonherited_model")
plt.savefig("im5_new.pdf")
plt.show()

data_final.plot(x='Population_size', y='Invasion_probability')
plt.ylabel('Invasion probabilities')
plt.title("Invasion probabilites distribution in beta_gom_time_dep_stochastic_nonherited_model")
plt.savefig("im6_new.pdf")
plt.show()
        
