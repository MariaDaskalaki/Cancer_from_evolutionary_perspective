#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 15:44:18 2020

@author: maria
"""
################################ Code for plotting the average competing cancerous cells after invasion across generations #######################

import csv
from itertools import chain
from itertools import accumulate
import re
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import numpy as np

path="/home/maria/CA_N500_gens5000_constant_gom_time_dep_stochastic_non_herited"
nruns=100
cancerous_fixation_list=[]
count_list=[]

path_for_population_size=path+"/results1/Number of mutants.txt"
with open (path_for_population_size, "r") as f_pop:
    data=f_pop.readlines()
    initial_population_size=int(re.sub("[^0-9]", "", data[0]))
    print("The population size is {}".format(initial_population_size))

    

for i in range(nruns):
    path_for_run_fixation=path+"/results{}/Number of mutants.txt".format(i+1)
    path_for_run_competition=path+"/results{}/Frequencies_CA_mutant.txt".format(i+1)
    
    with open (path_for_run_competition, "r") as main_file:
        main_data=main_file.readlines()
        generation_list=[int(line.split()[1]) for line in main_data]
        genot_freq_list=[line.split()[2:] for line in main_data]
        #print (genot_freq_list)
        count_list_run=[len(i) for i in genot_freq_list]
        count_list.append(count_list_run)
        
    with open (path_for_run_fixation, "r") as file:
        new_data=file.readlines()[1:]
        for line in new_data:
            find=False
            if int(line.split()[2])==initial_population_size:
                find=True
            if find==True:
                gen_invasion=int(line.split()[1])
                cancerous_fixation_list.append(gen_invasion)
                print ("The invasion has happened in generation {}".format(gen_invasion)) # find the generation where invasion of cancerous cells has hapenned
                break
        if find==False:
            print ("There has been no invasion of cancerous")
#print (generation_list)
print ("The len of count list is",len(count_list))
print (cancerous_fixation_list)
mean_cancerous_fixation_time=sum(cancerous_fixation_list)/len(cancerous_fixation_list)
print ("The mean fixation time of cancerous for population size {} is {}".format(initial_population_size,mean_cancerous_fixation_time))
rounded_mean_cancerous_fixation_time=round(mean_cancerous_fixation_time) # round the invasion time of mutants.
print ("So the rounded mean fixation time of cancerous for population size {} is {}".format(initial_population_size,rounded_mean_cancerous_fixation_time))

print ("The len of generation list is", len(generation_list))
new_list=[j[i] for i in generation_list for j in count_list]
split_length=len(count_list)
count_list_final = [new_list[x:x+split_length] for x in range(0, len(new_list), split_length)] # split list according to the generations of each run
final_dictionary=dict(zip(generation_list,count_list_final)) # dictionary with generations as keys and number of cancerous as values.

new_dictionary={i:(sum(final_dictionary[i])/len(final_dictionary[i])) for i in final_dictionary.keys() if i>= rounded_mean_cancerous_fixation_time} # dictionary with generations as keys and number of cancerous as values but after the time of invasion
print (new_dictionary)

plt.plot(list(new_dictionary.keys()),list(new_dictionary.values()))
plt.xlabel("Generations")
plt.ylabel("Average number of cancerous genotypes")
plt.title("Competition for fixation of cancerous genotypes")
plt.show() # competition plot after the time of invasion
        
        
