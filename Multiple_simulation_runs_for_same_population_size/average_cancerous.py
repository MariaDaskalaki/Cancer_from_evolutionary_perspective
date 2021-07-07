#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 13:47:57 2020

@author: maria
"""
###################### Code for plotting the average number of mutants form 100 simulation runs and the tumor size curve from 100 simulation runs of the same population size #############

import csv
from itertools import chain
from itertools import accumulate
import re
import pandas as pd
import matplotlib.pyplot as plt
import fileinput
from collections import defaultdict
import networkx as nx
import matplotlib.lines as mlines

nruns=100

path="/home/maria/Diplomatiki_parousiasi/no_nei_deterministic"
path_for_run=[]
path_for_population_size=path+"/results1/Number of mutants.txt"

with open (path_for_population_size, "r") as f_pop:
    data=f_pop.readlines()
    initial_population_size=int(re.sub("[^0-9]", "", data[0]))
    print("The population size is {}".format(initial_population_size))

string_to_exclude="The total number of initial healthy cells are {}".format(initial_population_size)


for i in range(nruns):
    path_for_run.append(path+"/results{}/Number of mutants.txt".format(i+1))
    
with open ("/home/maria/Desktop/res_file.txt", "w") as fout, fileinput.input(i for i in path_for_run) as fin:
    for line in fin:
        if line.strip("\n") != string_to_exclude:
            fout.write(line)
            
dictionary=defaultdict(list)


with open ("/home/maria/Desktop/res_file.txt", "r") as new_file:
    for line in new_file:
        key=int(line.split()[1])
        value=int(line.split()[2])
        dictionary[key].append(value)
        
#print (dictionary)

with open ("/home/maria/Desktop/average_cancerous.txt", "w+") as final_file:
    final_file.write("Generation\tAverage_number_of_cancerous\n")
    for (keys, values) in dictionary.items():
        #print (len(values)) # it is equal to the nruns
        average_number=sum(values)/len(values)
        #print("For generation {} the average number of mutants is {}".format(keys, average_number))
        final_file.write("{}\t{}\n".format(keys,average_number))
        
G=nx.Graph()
G.add_nodes_from(range(initial_population_size))


data_final=pd.read_csv("/home/maria/Desktop/average_cancerous.txt", sep='\t', lineterminator='\n', header=0)
length=(len(data_final["Average_number_of_cancerous"]))
#print (data_final["Average_number_of_cancerous"][0])
draw_list=[data_final["Average_number_of_cancerous"][0], data_final["Average_number_of_cancerous"][1], data_final["Average_number_of_cancerous"][int(length/2)],data_final["Average_number_of_cancerous"][length-1]] # time to plot the number of mutants 
print(draw_list)
red_square = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
                          markersize=10, label='Canceerous')
blue_square = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
                          markersize=10, label='Healthy')
#print(draw_list)
title_list=["Initial cells state for {} population".format(initial_population_size), "First generation cells state for {} population".format(initial_population_size),"Middle cells state for {} population".format(initial_population_size) ,"Final cells state for {} population".format(initial_population_size)]
print("The average numbers of cancerous for 0, first, middle and final generations are",draw_list)

cmap= {}
count=0
    
for i in draw_list:

    for node in G:

        if node < i:

            cmap[node]='red'
        else:

            cmap[node]='blue'
        #     print(node)
    
    plt.title(title_list[count])
    nx.draw_circular(G, with_labels=False, node_color=cmap.values(),node_size=50)
    #plt.title(title_list[count])
    count+=1
    #plt.legend("cancerous")
    #plt.gca().legend(('cancerous','healthy'))
    plt.legend(handles=[blue_square, red_square], loc='upper right')
    #plt.legend("healthy")

    plt.axis('equal')
    plt.show() # cirvular graph for cells, blue color for healthy cells and red color for mutants
    


data_final.plot(x='Generation', y='Average_number_of_cancerous')
plt.ylabel('Average number of cancerous cells')
plt.title('Average distribution of cancerous cells in time for {} population'.format(initial_population_size))
plt.show() # average tumor size curve 
