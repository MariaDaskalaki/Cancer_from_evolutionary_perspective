#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 17:36:09 2020

@author: maria
"""

import csv
from itertools import chain
from itertools import accumulate
import re
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import numpy as np

'''

with open ("/home/maria/Desktop/population_vs_fixation_time.txt", "w+") as final_file:
        final_file.write("Model\tPopulation_size\tMean_fixation_time\tProbability_of_fixation\tMean_invasion_fixation_time\tInvasion_probability\tMean_fitness\tMean_Fitness_variance\n")


nruns=100
models_dictionary={"beta_model":["/home/maria/CA_N100_gens5000_beta_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N300_gens5000_beta_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N500_gens5000_beta_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N800_gens5000_beta_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N1000_gens5000_beta_gom_time_dep_stochastic_non_herited"], "linear_model":["/home/maria/CA_N100_gens5000_linear_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N300_gens5000_linear_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N500_gens5000_linear_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N800_gens5000_linear_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N1000_gens5000_linear_gom_time_dep_stochastic_non_herited"], "constant_model":["/home/maria/CA_N100_gens5000_constant_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N300_gens5000_constant_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N500_gens5000_constant_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N800_gens5000_constant_gom_time_dep_stochastic_non_herited", "/home/maria/CA_N1000_gens5000_constant_gom_time_dep_stochastic_non_herited"]}
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
        #path_for_population_size="/home/maria/CA_N100_gens1000_linear_linear_time_ind_deterministic_inherited/results1/Number of mutants.txt"
        path_for_population_size= path+"/results1/Number of mutants.txt"

        with open (path_for_population_size, "r") as f_pop:
            data=f_pop.readlines()
            initial_population_size=int(re.sub("[^0-9]", "", data[0]))
            population_size_list.append(initial_population_size)
            print("The population size is {}".format(initial_population_size))


        for i in range(nruns):
            path_for_run= path+"/results{}/Frequencies_CA.txt".format(i+1)
            path_for_run_cancerous=path+"/results{}/Number of mutants.txt".format(i+1)
            path_for_run_competition=path+"/results{}/Frequencies_CA_mutant.txt".format(i+1)
            path_for_run_variance=path+"/results{}/Average_fitness.txt".format(i+1)
            #path_for_run="/home/maria/CA_N100_gens1000_linear_linear_time_ind_deterministic_inherited/results{}/Frequencies_CA.txt".format(i+1)
            generation_count=0
            with open (path_for_run, "r") as f:
                data=f.readlines()
                generation_list=[int(line.split()[1]) for line in data]
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
             
                for i in freq_list_new:
                    find=False
                    generation_count+=1
                    for j in i:
                        if j==1.0:
                            find=True
                            #print("There has been fixation in generation ",generation_count)
                    #break
                    if find==True:
                        #fixation_count+=1
                        print("There has been fixation in generation ",generation_count-1)
                        fixation_list.append(generation_count-1)
                        key=initial_population_size
                        value=generation_count-1
                        default_dictionary[key].append(value)
                        break
                    #else:
                     #   print("There has been no fixation")
                      #  break
                if find==False:
                    print("There has been no fixation")

            print ("---------------Start of cancerous study--------------")


            with open (path_for_run_cancerous, "r") as file:
                    main_data=file.readlines()[1:]
                    for line in main_data:
                        find=False
                        if int(line.split()[2]) == initial_population_size:
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
            total_average_means=(sub_size*sum(data_variance[1]))/total_population_size
            variance=1/total_population_size*(sum((data_variance[1]-total_average_means)**2 + data_variance[2]**2))
            #print (variance)
            variance_list.append(variance)
            means_fitness_list.append(total_average_means)

        mean_variance_accross_runs=sum(variance_list)/len(variance_list)
        mean_fitness_across_runs=sum(means_fitness_list)/len(means_fitness_list)
        print (cancerous_fixation_list)
        print ("There have been {} fixation of cancerous".format(len(cancerous_fixation_list)))
        mean_cancerous_fixation_time=sum(cancerous_fixation_list)/len(cancerous_fixation_list)
        print ("The mean fixation time of cancerous for population size {} is {}".format(initial_population_size,mean_cancerous_fixation_time))
        rounded_mean_cancerous_fixation_time=round(mean_cancerous_fixation_time)
        print ("So the rounded mean fixation time of cancerous for population size {} is {}".format(initial_population_size,rounded_mean_cancerous_fixation_time))
        invasion_probability=len(cancerous_fixation_list)/nruns
        print ("The invasion probability for population size {} is {}".format(initial_population_size, invasion_probability))
        print ("The len of count list is", len(count_list))
        #print ("The len of generation list is",len(generation_list))
        #print (count_list)
        new_list=[j[i] for i in generation_list for j in count_list]
        split_length=len(count_list)
        count_list_final = [new_list[x:x+split_length] for x in range(0, len(new_list), split_length)]
        final_dictionary=dict(zip(generation_list,count_list_final))
        new_dictionary={i:(sum(final_dictionary[i])/len(final_dictionary[i])) for i in final_dictionary.keys() if i>= rounded_mean_cancerous_fixation_time}
        competition_dictionaries_list.append(new_dictionary)

        print(fixation_list)
        #print("There are {} fixations".format(fixation_count))
        print("There are {} fixations".format(len(fixation_list)))

        mean_fixation_time=sum(fixation_list)/len(fixation_list)
        probability_of_fixation=len(fixation_list)/nruns
        print("The mean fixation time is {}".format(mean_fixation_time))
        print("The probability of fixation is {}".format(probability_of_fixation))



        with open ("/home/maria/Desktop/population_vs_fixation_time.txt", "a+") as final_file:
            #final_file.write("Population_size\tMean_fixation_time\n")
            final_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(model,initial_population_size, mean_fixation_time, probability_of_fixation, mean_cancerous_fixation_time, invasion_probability, mean_fitness_across_runs, mean_variance_accross_runs))

'''

df=pd.read_csv("/home/maria/Diplomatiki_parousiasi/population_vs_fixation_time.txt", sep='\t', lineterminator='\n', header=0)
new_df_100=df[df['Population_size']==100]
new_df_300=df[df['Population_size']==300]
new_df_500=df[df['Population_size']==500]
new_df_800=df[df['Population_size']==800]
new_df_1000=df[df['Population_size']==1000]
#print (new_df)

sns.lineplot(data=df, x='Population_size', y='Mean_fixation_time', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Generations")
plt.title("Mean fixation time distribution")
plt.show()

sns.lineplot(data=df, x='Population_size', y='Probability_of_fixation', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Fixation probabilities")
plt.title("Fixation probabilities distribution")
plt.show()

sns.lineplot(data=df, x='Population_size', y='Mean_invasion_fixation_time', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Generations")
plt.title("Mean invasion time of cancerous cells distribution")
plt.show()

sns.lineplot(data=df, x='Population_size', y='Invasion_probability', hue='Model')
plt.xlabel("Population size")
plt.ylabel("Invasion probabilities")
plt.title("Invasion probabilities of cancerous cells distribution")
plt.show()

sns.lineplot(data=df, x='Mean_Fitness_variance', y='Mean_fixation_time', hue='Model')
plt.xlabel("Mean fitness variances")
plt.ylabel("Mean fixation times")
plt.title("Variance vs Fixation time in different models")
plt.show()

sns.lineplot(data=new_df_100, x='Mean_Fitness_variance', y='Mean_fixation_time', label="100pop")
plt.xlabel("Mean fitness variances")
plt.ylabel("Mean fixation times")
plt.title("Variance vs Fixation time in different models for 100 population size")
plt.show()

sns.lineplot(data=new_df_300, x='Mean_Fitness_variance', y='Mean_fixation_time', label="300pop")
#sns.lineplot(data=new_df_300, x='Mean_Fitness_variance', y='Mean_fixation_time')
plt.xlabel("Mean fitness variances")
plt.ylabel("Mean fixation times")
plt.title("Variance vs Fixation time in different models for 300 population size")
plt.show()

sns.lineplot(data=new_df_500, x='Mean_Fitness_variance', y='Mean_fixation_time', label="500pop")
plt.xlabel("Mean fitness variances")
plt.ylabel("Mean fixation times")
plt.title("Variance vs Fixation time in different models for 100 population size for 500 population size")
plt.show()

sns.lineplot(data=new_df_800, x='Mean_Fitness_variance', y='Mean_fixation_time', label="800pop")
plt.xlabel("Mean fitness variances")
plt.ylabel("Mean fixation times")
plt.title("Variance vs Fixation time in different models for 100 population size for 800 population size")
plt.show()

sns.lineplot(data=new_df_1000, x='Mean_Fitness_variance', y='Mean_fixation_time', label="1000pop")
plt.xlabel("Mean fitness variances")
plt.ylabel("Mean fixation times")
plt.title("Variance vs Fixation time in different models for 100 population size for 1000 population size")
plt.show()




    
