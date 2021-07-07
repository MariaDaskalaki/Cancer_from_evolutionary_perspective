#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 03:01:02 2020

@author: maria
"""

################### Code for plotting discrete distributions of frequencies for each mutation of each genome position #####################

import csv
from itertools import chain
from itertools import accumulate
import re
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#import numpy as np
import time


nruns=100
dfs_lista=[]

generations_lista=[]
for i in range(1,nruns+1):
    path="/home/maria/Diplomatiki_parousiasi/parousiasi_deterministic_new/results{}/Frequencies_mutations_CA.txt".format(i)
    df=pd.read_csv(path,sep="\t",header=None)
    generations_lista.append(len(df)-1)
    #if i==2:
     #   break
print (generations_lista)
max_generation=max(generations_lista)
print (max_generation) # the maximum generation of fixation from 100 simulation runs

for i in range(1,nruns+1):
    path="/home/maria/Diplomatiki_parousiasi/parousiasi_deterministic_new/results{}/Frequencies_mutations_CA.txt".format(i)
    df=pd.read_csv(path,sep="\t",header=None)
    df.drop(df.columns[-1], axis=1, inplace=True)
    df[0]=df[0].str.split().str.get(1).astype(float)
    selected_df=df
    #print (selected_df)
    transposed_df=selected_df.T
    #print (transposed_df)
    transposed_df.columns = transposed_df.iloc[0]
    transposed_df = transposed_df.reindex(transposed_df.index.drop(0)).reset_index(drop=True)
    transposed_df.columns.name = None
    transposed_df = transposed_df.apply(lambda x: x.str.split(r'\:').str.get(1).astype(float), axis=1)
    #print (transposed_df)
    new_column=pd.Series(transposed_df[transposed_df.columns[-1]])
    if max_generation!=len(selected_df)-1:
        extra_df=pd.concat([new_column]*(max_generation-len(selected_df)+1), axis=1) # add columns, until the max generation, with the same frequencies as the last column of each run
        new_big_df=pd.concat([transposed_df, extra_df],ignore_index=True, axis=1)
    else:
        new_big_df=transposed_df
    #print (new_big_df)
    dfs_lista.append(new_big_df)
    #if i==2:
     #   break
averages = pd.concat([each.stack() for each in dfs_lista],axis=1).apply(lambda x:x.mean(),axis=1).unstack() # find the average frequency of each mutation across generations
print (averages)
#print (averages[111])
    
    
#averages = pd.concat([each.stack() for each in dfs_lista],axis=1).apply(lambda x:x.mean(),axis=1).unstack()
mut_dict_final=averages.to_dict('dict')
#print (mut_dict_final)
print ("The averages freqs are")
print (averages)
gen_number=len(averages.columns)
print ("The number of generations is",gen_number)
wanted_gen_list=[float(i) for i in range(0,len(averages.columns), 300)]
print (wanted_gen_list)
subset_mut_dict_final=dict((k,mut_dict_final[k]) for k in wanted_gen_list)
new_dict={}

for k,v in subset_mut_dict_final.items():
    new_dict[k]=list(v.values())
    
subset_df=pd.DataFrame([(key, var) for (key, L) in new_dict.items() for var in L], columns=['Generations', 'Frequencies'])
subset_df['Generations']=subset_df['Generations'].astype(int)
subset_df['Coordinates'] = list(zip(subset_df.Generations, subset_df.Frequencies)) # dataframe with average frequency of each mutation in each generation

z=averages.values
#print (averages.values)

length=len(averages.columns)

fig, ax = plt.subplots(figsize=(6,6))
plt.imshow(z,origin='lower', aspect='auto',extent=[0,length, 0, len(averages.index)])
plt.colorbar()
plt.title("Mutation Frequencies heatmap density plot")
plt.xlabel("Generations")
plt.ylabel("Genome Positions")

averages=averages.transpose()
#print(df)

averages.plot(figsize=(10,6), xticks=range(0,len(averages),300)).legend(title='', bbox_to_anchor=(1, 1))
plt.title("All mutations evolution")
plt.xlabel("Generations")
plt.ylabel("Frequencies")
plt.show() 

fig, ax = plt.subplots(figsize=(14,9))
g = sns.violinplot('Frequencies','Generations', data=subset_df, linewidth=0.2, orient="h", cut=0)
#g.set_xscale('log', basex=10)
ax.set(ylabel='Generations', xlabel='Frequencies')
plt.title("Violin density plot for mutations positions freqs")
plt.show()

subset_df.plot(x ='Generations', y='Frequencies', kind = 'scatter')
plt.title("Discrete distribution of mutations positions freqs accross gens")
plt.show()

ax2 = sns.swarmplot(x="Generations", y="Frequencies", data=subset_df)
plt.title("Colored scatterplot of discrete distribution of mutations positions freqs accross gens")
plt.show()



