#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 17:36:09 2020

@author: maria
"""
############### Code for plotting mean fixation and inivasion time and probability of fixation and invasion for different fitness site scenarios #############################

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
        final_file.write("Model\tPopulation_size\tMean_fixation_time\tProbability_of_fixation\tMean_invasion_fixation_time\tInvasion_probability\tMean_fitness\tMean_Fitness_variance\n") # create file with the statistics of all simulation runs for each fitness scite scenario for all the population sizes


nruns=100
models_dictionary={"beta_model":["/home/daskalaki/synology/daskalaki/100pop_new", "/home/daskalaki/synology/daskalaki/300pop_new", "/home/daskalaki/synology/daskalaki/500pop_new", "/home/daskalaki/synology/daskalaki/800pop_new", "/home/daskalaki/synology/daskalaki/1000pop_new"], "linear_model":["/home/daskalaki/synology/daskalaki/100pop_linear_new", "/home/daskalaki/synology/daskalaki/300pop_linear_new", "/home/daskalaki/synology/daskalaki/500pop_linear_new", "/home/daskalaki/synology/daskalaki/800pop_linear_new", "/home/daskalaki/synology/daskalaki/1000pop_linear_new"], "constant_model":["/home/daskalaki/synology/daskalaki/100pop_con_new", "/home/daskalaki/synology/daskalaki/300pop_con_new", "/home/daskalaki/synology/daskalaki/500pop_con_new", "/home/daskalaki/synology/daskalaki/800pop_con_new", "/home/daskalaki/synology/daskalaki/1000pop_con_new"],"sin_model":["/home/daskalaki/synology/daskalaki/100pop_alter_new", "/home/daskalaki/synology/daskalaki/300pop_alter_new", "/home/daskalaki/synology/daskalaki/500pop_alter_new", "/home/daskalaki/synology/daskalaki/800pop_alter_new", "/home/daskalaki/synology/daskalaki/1000pop_alter_new"]}
#path_for_population_size="/home/maria/CA_N100_gens1000_linear_linear_time_ind_deterministic_inherited_test/results1/Number of mutants.txt"
#fixation_count=0
#fixation_list=[]
for (models,folders) in models_dictionary.items():


    path_list=folders
    default_dictionary=defaultdict(list)
    default_cancerous_dictionary=defaultdict(list)
    population_size_list=[]
    competition_dictionaries_list=[]
    model=models

    for path in path_list:
        variance_list=[]
        means_fitness_list=[]
        fixation_list=[]
        cancerous_fixation_list=[]
        count_list=[]
        path_for_population_size= path+"/results1/Number of mutants.txt"

        with open (path_for_population_size, "r") as f_pop:
            data=f_pop.readlines()
            initial_population_size=int(re.sub("[^0-9]", "", data[0]))
            population_size_list.append(initial_population_size)
            print("The population size is {}".format(initial_population_size)) #initial population


        for i in range(nruns):
            path_for_run= path+"/results{}/Frequencies_CA.txt".format(i+1)
            path_for_run_cancerous=path+"/results{}/Number of mutants.txt".format(i+1)
            path_for_run_competition=path+"/results{}/Frequencies_CA_mutant.txt".format(i+1)
            path_for_run_variance=path+"/results{}/Average_fitness.txt".format(i+1)
            generation_count=0
            with open (path_for_run, "r") as f:
                data=f.readlines()
                generation_list=[int(line.split()[1]) for line in data]
                fixation_generation=generation_list[-1] # time of fixation for each run
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
                #print(freq_list_new)
                fixation_list.append(fixation_generation)
                key=initial_population_size
                value=fixation_generation
                default_dictionary[key].append(value)
              



            with open (path_for_run_cancerous, "r") as file:
                    main_data=file.readlines()[1:]
                    for line in main_data:
                        find=False
                        if int(line.split()[2]) == initial_population_size:
                            find=True
                        if find==True:
                            fixed_generation=int(line.split()[1])
                            cancerous_fixation_list.append(fixed_generation) # time of invasion of mutants for each run
                            key_cancerous=initial_population_size
                            value_cancerous=fixed_generation
                            default_cancerous_dictionary[key_cancerous].append(value_cancerous)
                            print("The fixation of cancerous has happened in generation {}".format(fixed_generation))
                            break
                    if find==False:
                        print ("The has been no fixation of cancerous")

            with open (path_for_run_competition, "r") as main_file:
                main_data=main_file.readlines()
                generation_list=[int(line.split()[1]) for line in main_data]
                genot_freq_list=[line.split()[2:] for line in main_data]
                #print (genot_freq_list)
                count_list_run=[len(i) for i in genot_freq_list]
                #print ("The count_list is", count_list_run)
                count_list.append(count_list_run) #list of lists with number of competitive cells for each generation of each run
                
            data_variance=pd.read_csv(path_for_run_variance, sep='\t', lineterminator='\n', header=None)
            generations=len(data_variance)
            generations_list=[j for j in range(generations)]
            sub_size=initial_population_size
            total_population_size=generations*sub_size
            total_average_means=(sub_size*sum(data_variance[1]))/total_population_size # fitness for each run 
            variance=1/total_population_size*(sum((data_variance[1]-total_average_means)**2 + data_variance[2]**2)) # fitness variance for each run
            #print (variance)
            variance_list.append(variance) #list of fitness variances for each run
            means_fitness_list.append(total_average_means) # list of mean fitness for each run

        mean_variance_accross_runs=sum(variance_list)/len(variance_list) # mean fitness variace from all the simulation runs of the same fitness site model of each population
        mean_fitness_across_runs=sum(means_fitness_list)/len(means_fitness_list) # mean fitness from all the simulation runs of the same fitness site model of each population
        print (cancerous_fixation_list)
        print ("There have been {} fixation of cancerous".format(len(cancerous_fixation_list)))
        mean_cancerous_fixation_time=sum(cancerous_fixation_list)/len(cancerous_fixation_list) # mean invasion time of the same fitness site model for each population
        print ("The mean fixation time of cancerous for population size {} is {}".format(initial_population_size,mean_cancerous_fixation_time))
        rounded_mean_cancerous_fixation_time=round(mean_cancerous_fixation_time) # rounded mean invasion time
        print ("So the rounded mean fixation time of cancerous for population size {} is {}".format(initial_population_size,rounded_mean_cancerous_fixation_time))
        invasion_probability=len(cancerous_fixation_list)/nruns # invasion probability 
        print ("The invasion probability for population size {} is {}".format(initial_population_size, invasion_probability))
        print ("The len of count list is", len(count_list))
        split_length=len(count_list)
        print(fixation_list)
        print("There are {} fixations".format(len(fixation_list)))

        mean_fixation_time=sum(fixation_list)/len(fixation_list) # mean fixation time of the same fitenss site model for each population 
        probability_of_fixation=len(fixation_list)/nruns # fixation probability 
        print("The mean fixation time is {}".format(mean_fixation_time))
        print("The probability of fixation is {}".format(probability_of_fixation))



        with open ("/home/daskalaki/synology/daskalaki/population_vs_fixation_time.txt", "a+") as final_file:
            #final_file.write("Population_size\tMean_fixation_time\n")
            final_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(model,initial_population_size, mean_fixation_time, probability_of_fixation, mean_cancerous_fixation_time, invasion_probability, mean_fitness_across_runs, mean_variance_accross_runs))
            

df=pd.read_csv("/home/daskalaki/synology/daskalaki/population_vs_fixation_time.txt", sep='\t', lineterminator='\n', header=0)

sns.lineplot(data=df, x='Population_size', y='Mean_fixation_time', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Generations")
plt.title("Mean fixation time distribution")
plt.savefig("/home/daskalaki/synology/daskalaki/Images_new_new/Mean_fixation_time_distribution.png")
plt.show()

sns.lineplot(data=df, x='Population_size', y='Probability_of_fixation', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Fixation probabilities")
plt.title("Fixation probabilities distribution")
plt.savefig("/home/daskalaki/synology/daskalaki/Images_new_new/Fixation_probabilities_distribution.png")
plt.show()

sns.lineplot(data=df, x='Population_size', y='Mean_invasion_fixation_time', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Generations")
plt.title("Mean invasion time of cancerous cells distribution")
plt.savefig("/home/daskalaki/synology/daskalaki/Images_new_new/Mean_invasion_time_of_cancerous_cells_distribution.png")
plt.show()

sns.lineplot(data=df, x='Population_size', y='Invasion_probability', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Invasion probabilities")
plt.title("Invasion probabilities of cancerous cells distribution")
plt.savefig("/home/daskalaki/synology/daskalaki/Images_new_new/Invasion_probabilities_of_cancerous_cells_distribution.png")
plt.show()

sns.lineplot(data=df, x='Mean_Fitness_variance', y='Mean_fixation_time', hue='Model')
plt.xlabel("Mean fitness variances")
plt.ylabel("Mean fixation times")
plt.title("Variance vs Fixation time in different models")
plt.savefig("/home/daskalaki/synology/daskalaki/Images_new_new/Variance_vs_Fixation_time_in_different_models.png")
plt.show()
    
