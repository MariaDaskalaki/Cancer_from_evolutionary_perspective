#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 01:23:44 2020

@author: maria
"""

###################### Code for plotting the mutants and the tumor size curve over time ###############################

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import pandas as pd
import re
import networkx as nx
import matplotlib.lines as mlines


with open ("Number of mutants.txt","r") as f:
    data=f.readlines()
    mutants_data=data[1:]
    initial_population_size=int(re.sub("[^0-9]", "", data[0]))
    #print(initial_population_size)
    x=[line.split()[0]for line in mutants_data] #generation
    y=[line.split()[1]for line in mutants_data] #number of generations
    gen_list=[int(i) for i in y]
    #print(len(gen_list))
    mutants_number=[int(line.split()[2]) for line in mutants_data]
    final_dict=dict(zip(gen_list,mutants_number))
    #print(final_dict)
    
    draw_list=[mutants_number[0], mutants_number[int(len(gen_list)/len(gen_list))], mutants_number[int(len(gen_list)/(int(len(gen_list)/12)))], mutants_number[int(len(gen_list)/2)], mutants_number[-1]] # times to plot the mutants of the circular graph
    title_list=["Initial cells state", "First generation cell state", "Before Middle cells state", "Middle cells state",  "Final cells state"]
    red_square = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
                          markersize=10, label='Canceerous')
    blue_square = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
                          markersize=10, label='Healthy')
    #print(draw_list)


    G=nx.Graph()
    G.add_nodes_from(range(initial_population_size))

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
        nx.draw_circular(G, with_labels=False, node_color=cmap.values(),node_size=50) # plot the circular graph of cells, blue color for healthy and red color for mutants
        #plt.title(title_list[count])
        count+=1
        #plt.legend("cancerous") plt.title(title_list[count])
        #plt.gca().legend(('cancerous','healthy'))
        plt.legend(handles=[blue_square, red_square], loc='upper right')
        #plt.legend("healthy")

        plt.axis('equal')
        plt.show()
        
    plt.plot(gen_list,mutants_number)
    #plt.xticks(np.arange(0,len(gen_list),100))
    plt.xlabel("Generations")
    plt.ylabel("Number of mutants")
    plt.title("Tumor size")
    plt.show() # tumor size curve plot
    
