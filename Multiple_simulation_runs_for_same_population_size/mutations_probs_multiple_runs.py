#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 23:15:47 2020

@author: maria
"""
#################### Code for plotting the probability of fixation for each genome position from 100 simulation runs of the same population size ################

import itertools
import matplotlib.pyplot as plt
import pandas as pd

nruns=100
fixed_positions_that_matter=[]
fixed_positions=[]
fixed_positions_matter_probs={}
fixed_positions_probs={}

for i in range(nruns):
    path="/home/maria/Diplomatiki_parousiasi/parousiasi_deterministic_new/results{}/Positions_that_matter.txt".format(i+1)
    
    
    with open(path,"r") as f:
        data = f.readlines()
        tail = data[-2:]
        #print (tail)
        a= tail[0].split(":")[1].split("\t")[:-1]
        b=tail[1].split(":")[1].split("\t")[:-1]
        #print (b)
        fixed_positions_that_matter.append(a)
        fixed_positions.append(b) # list with the fixed positions of each run
#print (fixed_positions_that_matter)

flatten_fixed_positions_that_matter = list(itertools.chain(*fixed_positions_that_matter))
flatten_fixed_positions=list(itertools.chain(*fixed_positions)) # returns a flattened list with all the fixed positions across 100 simulation runs
#print (flatten_fixed_positions_that_matter)

for i in flatten_fixed_positions_that_matter:
    fixed_positions_matter_probs[int(i)]=flatten_fixed_positions_that_matter.count(i)/nruns #probability of fixation for the positions with significance
    
for j in flatten_fixed_positions:
    fixed_positions_probs[int(j)]=flatten_fixed_positions.count(j)/nruns # probability of fixation for each genome position

fixed_positions_matter_probs=dict(sorted(fixed_positions_matter_probs.items()))
fixed_positions_probs=dict(sorted(fixed_positions_probs.items()))
print ("The probs for fixed positions that matter are")
print (fixed_positions_matter_probs)
print ("The probs for all the fixed positions are")
print (fixed_positions_probs)

colors = []
for key in fixed_positions_probs.keys(): # keys are the names of the boys
    if key in range(0,15):
        colors.append('r')
    else:
        colors.append('b')

#bar(ind,num,width,color=colors)

plt.bar(fixed_positions_matter_probs.keys(), fixed_positions_matter_probs.values())
plt.xticks(list(fixed_positions_matter_probs.keys()))
plt.xlabel("Genome positions")
plt.ylabel("Probabilities")
plt.title("Fixation probabilities distribution of genome positions mutations from those that matter")
plt.show()

plt.bar(fixed_positions_probs.keys(), fixed_positions_probs.values(),color=colors, label="Positions that matter")
plt.xticks(list(fixed_positions_probs.keys()))
plt.xlabel("Genome positions")
plt.ylabel("Probabilities")
plt.legend()
plt.title("Fixation probabilities distribution of genome positions mutations from those that matter and not")
plt.show()

data_weights=pd.read_csv("/home/maria/Diplomatiki_parousiasi/parousiasi_deterministic_new/results1/Position_weights.txt", sep="\t", header=None)
data_weights.columns=["Positions", "Weights"]
print (data_weights)

data_weights.plot(x="Positions", y="Weights", kind='scatter')
plt.xticks(list(data_weights["Positions"]))
plt.title("Positions weights")
plt.show() # plot of weights for from gaussian distribution for each genome position
