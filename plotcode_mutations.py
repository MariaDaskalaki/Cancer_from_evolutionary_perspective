#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:44:15 2020

@author: maria
"""

import matplotlib.pyplot as plt
from itertools import chain
from itertools import accumulate
from collections import defaultdict
import pandas as pd
import itertools

with open("Frequencies_mutations.txt",'r') as f:
    data=f.readlines()
    x=[line.split()[0]for line in data] #generation
    y=[line.split()[1]for line in data] #number of generations
    mut_freq=[line.split()[2:] for line in data] # lists of each generation with mutation:frequencies
    gen_list=["{}".format(i) for i in range(len(mut_freq))]
    gen_list=[float(i) for i in gen_list]
    mut_freq_lista= list(chain.from_iterable(mut_freq)) # whole list of mutations:frequencies
    mut=[i.split(":")[0] for i in mut_freq_lista] # list with all the mutations
    mut_enum=[i[0] for i in enumerate(mut)] # enumerated mutations
    freq=[i.split(":")[1] for i in mut_freq_lista] #list with all str(frequencies)
    float_freq=[float(i) for i in freq] #list with all float(frequencies)
    string_floats_freq=[str(i) for i in float_freq] #str(float)frequencies
    lista=[i+":"+j for i,j in zip(mut,string_floats_freq)] # list of lists for each generation, strings genotypes:str(float) frequencies
    split_length=[len(i) for i in mut_freq] #list of how many elements in order to return to the initial list of lists
    mut_freq_new=[lista[x - y: x] for x, y in zip(accumulate(split_length), split_length)] #return to the initial genot_freq_list but with str(floats) frequencies
    mut_freq_split=[j.strip().split(":") for i in mut_freq_new for j in i] # multiple lists of mutations,freqs without the :
    whole_mut_freq_split=list(chain.from_iterable(mut_freq_split)) # one list of mutations,freqs without the :
    new_split_length=[2*len(i) for i in mut_freq] 
    mut_freq_split_final=[whole_mut_freq_split[x - y: x] for x, y in zip(accumulate(new_split_length), new_split_length)] # break to the initial number of lists, one list for each generation mut,freq
    mut_freq_split_final_iteration=[iter(i) for i in mut_freq_split_final]
    mut_freq_split_final_tuples=[tuple(zip(i,i)) for i in mut_freq_split_final_iteration]
    mut_freq_split_tuples_list=[i for i in mut_freq_split_final_tuples] # tuples for each generation, so 3 master tuples
    dictionaries_mut_freq_list=[] # list of three dictionaries (one for each generation) with mut as key an frequency as value !!! not duplicate keys for each generation
    for i in mut_freq_split_tuples_list:
        dictt=defaultdict(list)
        #print(i)
        for k,v in i:
            dictt[k].append(v)
            #print(dictt[k]) εάν το βάλω εδώ και το αρχικοποιήσω και στην αρχή τοτε εκτυπώνει ένα λεξικο με mutations as keys και μία λίστα απ όλα τα dreqs του συγκεκριμένου genotype για όλες τις γενιές. 
        dictionaries_mut_freq_list.append(dict(dictt))
    dictionaries_mut_freq_list_new=[] #list with three dictionaries, one for each generation,with genotypes as keys an float frequencies as values
    for i in dictionaries_mut_freq_list:
        dict1={}
        for k,v in i.items():
            for j in v:
                dict1[k]=float(j)
        dictionaries_mut_freq_list_new.append(dict1)
    mut_dict_new=dict(zip(gen_list,dictionaries_mut_freq_list))# full dictionary with generations as keys and for each generation key a dictionary as values with mutations:frequencies with not duplicate keys as mutations.
    mut_dict_final=dict(zip(gen_list,dictionaries_mut_freq_list_new)) #full dictionary with generations as keys and dictionaries of mutations and float frequencies as values.
    mut_dict=dict(zip(gen_list,(mut_freq_new))) # dictionary generation i: list of genotypes:frequencies
    half_freq_thres=(max(float_freq))/2
    freq_thres=(max(float_freq))*0.9
    enum_mut_freq_dict=dict(zip(mut_enum,float_freq)) # dictionary enumerated genotype:float frequency
    mut_dict=dict(zip(mut_enum,mut)) # dictionary enumerated mutation:str(mutation)
    matter_muts=[mut_dict[k] for k,v in enum_mut_freq_dict.items() if v>freq_thres]
    matter_muts_unique=list(set(matter_muts))
    #print (mut_dict_final)
    print(matter_muts_unique)
   
    
   
    
    
    df=pd.DataFrame(mut_dict_final)
    df=df.transpose()
    #print(df)
    
    df.plot(figsize=(10,6), xticks=range(0,len(gen_list),1000)).legend(title='', bbox_to_anchor=(1, 1))
    plt.title("All mutations evolution")
    plt.xlabel("Generations")
    plt.ylabel("Frequencies")
    plt.show() 
    
    df_new=df.loc[:,('{}'.format(i) for i in matter_muts_unique)]
    #print(df_new)
    df_new.plot(figsize=(10,6), xticks=range(0, len(gen_list),1000)).legend(title='', bbox_to_anchor=(1, 1))
    #df_new.plot(figsize=(6000,6), xticks=range(0,len(gen_list),1000)).legend(title='', bbox_to_anchor=(1, 1))
    plt.title("High frequency mutations evolution")
    plt.xlabel("Generations")
    plt.ylabel("Frequencies")
    plt.show()
    