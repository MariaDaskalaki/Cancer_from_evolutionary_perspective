#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 12:44:23 2020

@author: maria
"""
################3 Code for creating the statistical_analysis.txt and statistical_analysis_anova.txt with the fixation times of each one of the 100 simulation runs of each fitness site model #############

import csv
from itertools import chain
from itertools import accumulate
import re
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import numpy as np
import csv
from csv import DictWriter
from itertools import chain, repeat

nruns=100
path_list=["/home/daskalaki/synology/daskalaki/100pop_linear_new", "/home/daskalaki/synology/daskalaki/100pop_new", "/home/daskalaki/synology/daskalaki/100pop_con_new", "/home/daskalaki/synology/daskalaki/100pop_alter_new"]
#default_dictionary=defaultdict(list)
models_list=["linear_model", "beta_model", "constant_model", "sin_model"]
models_times_dictionary={}
index_models=0


for path in path_list:
    fixation_list=[]
    path_for_population_size= path+"/results1/Number of mutants.txt"
    
    with open (path_for_population_size, "r") as f_pop:
        data=f_pop.readlines()
        initial_population_size=int(re.sub("[^0-9]", "", data[0]))
        #population_size_list.append(initial_population_size)
        print("The population size is {}".format(initial_population_size))
    
    for i in range(nruns):
        path_for_run=path+"/results{}/Frequencies_CA.txt".format(i+1)
        generation_count=0
        with open (path_for_run, "r") as f:
            data=f.readlines()
            generation_list=[int(line.split()[1]) for line in data]
            fixation_generation=generation_list[-1] # generation where fixation happened 
            #print (generation_list)
            genot_freq_list=[line.split()[2:] for line in data]
            #print(genot_freq_list)
            #print("There has been fixation in generation ",count)
            genot_freq_list_final=list(chain.from_iterable(genot_freq_list))
            #print(genot_freq_list_final)
            genot_list=[i.split(':')[0] for i in genot_freq_list_final]
            #print (genot_list)
            freq_list=[float(i.split(':')[1]) for i in genot_freq_list_final]
            #print(freq_list)
            split_length=[len(i) for i in genot_freq_list]
            #print(split_length)
            freq_list_new=[freq_list[x - y: x] for x, y in zip(accumulate(split_length), split_length)]
            fixation_list.append(fixation_generation) # list of fixation times from each one of the 100 simulation runs of each fitenss site scenario
            #print(freq_list_new)
            
    #print (fixation_list)
    models_times_dictionary[models_list[index_models]]=fixation_list # dictionary with fitenss site scenario as keys and list of fixation times from 100 simulation runs of the corresponding scenrio as values
    index_models+=1
print (models_times_dictionary)



df = pd.DataFrame({ key:pd.Series(value) for key, value in models_times_dictionary.items() })
df.to_csv('/home/daskalaki/synology/daskalaki/statistical_analysis.txt', sep='\t',index=False)



data_anova = pd.DataFrame.from_dict(models_times_dictionary, orient='index')
data_anova = data_anova.stack().to_frame().reset_index().drop('level_1', axis=1)
data_anova.columns = ['Model', 'Time']
print (data_anova)

data_anova.to_csv("/home/daskalaki/synology/daskalaki/statistical_analysis_anova.txt", sep='\t', index=False)
