#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 20:22:23 2020

@author: maria
"""

############################## Code for plotting the average number of competitve cancerous cells for differnt population sizes across generations ##################

import pandas as pd
import matplotlib.pyplot as plt
import fileinput
from collections import defaultdict
#import networkx as nx
import re

nruns=100

path_list=["/home/daskalaki/synology/daskalaki/100pop_new","/home/daskalaki/synology/daskalaki/300pop_new","/home/daskalaki/synology/daskalaki/500pop_new","/home/daskalaki/synology/daskalaki/800pop_new","/home/daskalaki/synology/daskalaki/1000pop_new"]
whole_dictionary_list=[]
population_size_list=[]
#path="/home/maria/CA_N100_gens5000_linear_gom_time_dep_stochastic_non_herited"

for path in path_list:
    path_for_run=[]
    #dictionary={}
    path_for_population_size=path+"/results1/Number of mutants.txt"

    with open (path_for_population_size, "r") as f_pop:
        data=f_pop.readlines()
        initial_population_size=int(re.sub("[^0-9]", "", data[0]))
        population_size_list.append(initial_population_size) # initial population of helathy cells
        print("The population size is {}".format(initial_population_size))

    string_to_exclude="The total number of initial healthy cells are {}".format(initial_population_size)


    for i in range(nruns):
        path_for_run.append(path+"/results{}/Number of mutants.txt".format(i+1))

    with open ("/home/daskalaki/synology/daskalaki/res_file.txt", "w") as fout, fileinput.input(i for i in path_for_run) as fin: # res file for number of mutants of each generation for each run of each population size
        for line in fin:
            if line.strip("\n") != string_to_exclude:
                fout.write(line)

    dictionary=defaultdict(list)


    with open ("/home/daskalaki/synology/daskalaki/res_file.txt", "r") as new_file:
        for line in new_file:
            key=int(line.split()[1])
            value=int(line.split()[2])
            dictionary[key].append(value) # converts the res file into a dictionary with generations as keys and number of mutants as values

    #print (dictionary)

    with open ("/home/daskalaki/synology/daskalaki/average_cancerous_{}pop.txt".format(initial_population_size), "w+") as final_file:
        final_file.write("Generation\tAverage_number_of_cancerous\n")
        for (keys, values) in dictionary.items():
            #print (len(values)) # it is equal to the nruns
            average_number=sum(values)/len(values) # average number of mutants for each generation, from 100 simulation runs of the same population size
            #print("For generation {} the average number of mutants is {}".format(keys, average_number))
            final_file.write("{}\t{}\n".format(keys,average_number))

   

    data_final=pd.read_csv("/home/daskalaki/synology/daskalaki/average_cancerous_{}pop.txt".format(initial_population_size), sep='\t', lineterminator='\n', header=0)
    dictionary=dict(zip(data_final.Generation,data_final.Average_number_of_cancerous))
    whole_dictionary_list.append(dictionary) # list of dictionaries with generations as keys and average number of mutants as values 
    
index_count=0

for i in whole_dictionary_list:
        plt.plot(list(i.keys()),list(i.values()),label='population_size {}'.format(population_size_list[index_count]))
        index_count+=1
plt.legend(loc='upper right')
plt.xlabel("Generations")
plt.ylabel("Average number of cancerous cells")
plt.title("Average distribution of cancerous cells in time")
plt.savefig("im_average_cancerous_new.pdf")
plt.show() # plot the list of dictionaries


