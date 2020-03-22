#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:00:46 2020

@author: maria
"""

import matplotlib.pyplot as plt
from itertools import chain
from itertools import accumulate
from collections import defaultdict
import pandas as pd
import itertools




with open("Frequencies.txt",'r') as f:
    data=f.readlines()
    x=[line.split()[0]for line in data] #generation
    y=[line.split()[1]for line in data] #number of generations
    genot_freq=[line.split()[2:] for line in data] # lists of each generation with genotype:frequencies
    #print(genot_freq)
    #gen_list=["Generation {}".format(i) for i in range(len(genot_freq))] # list of strings generations
    gen_list=["{}".format(i) for i in range(len(genot_freq))]
    gen_list=[float(i) for i in gen_list]
    #print(gen_list)
    #gen_list_enum=[i for i in enumerate(gen_list)]
    genot_freq_lista= list(chain.from_iterable(genot_freq)) # all genotype:frequencies
    #print(genot_freq)
    genot=[i.split(":")[0] for i in genot_freq_lista] # all genotypes
    #print(genot)
    genot_enum=[i[0] for i in enumerate(genot)] # enumerated genotypes
    freq=[i.split(":")[1] for i in genot_freq_lista] # all frequencies
    float_freq=[float(i) for i in freq] # float all frequencies
    string_floats_freq=[str(i) for i in float_freq] #str(float)frequencies
    lista=[i+":"+j for i,j in zip(genot,string_floats_freq)] # list of lists for each generation, strings genotypes:str(float) frequencies
    #print(lista)
    split_length=[len(i) for i in genot_freq] #list of how many elements in order to return to the initial list of lists
    genot_freq_new=[lista[x - y: x] for x, y in zip(accumulate(split_length), split_length)] #return to the initial genot_freq_list but with str(floats) frequencies
    #print(genot_freq_new)
    genot_freq_split=[j.strip().split(":") for i in genot_freq_new for j in i] # multiple lists of genots,freqs without the :
    #print(genot_freq_split)
    whole_genot_freq_split=list(chain.from_iterable(genot_freq_split)) # one list of genots,freqs without the : 
    #print(whole_genot_freq_split)
    new_split_length=[2*len(i) for i in genot_freq] 
    genot_freq_split_final=[whole_genot_freq_split[x - y: x] for x, y in zip(accumulate(new_split_length), new_split_length)] # break to the initial number of lists, one list for each generation genot,freq
    #print(genot_freq_split_final)
    genot_freq_split_final_iteration=[iter(i) for i in genot_freq_split_final]
    #print(genot_freq_split_final_iteration)
    genot_freq_split_final_tuples=[tuple(zip(i,i)) for i in genot_freq_split_final_iteration]
    #print(genot_freq_slpit_final_tuples)
    genot_freq_split_tuples_list=[i for i in genot_freq_split_final_tuples] # tuples for each generation, so 3 master tuples
    #print(genot_freq_split_tuples_list)
    dictionaries_genot_freq_list=[] # list of three dictionaries (one for each generation) with genot as key an frequency as value !!! not duplicate keys for each generation
    for i in genot_freq_split_tuples_list:
        dictt=defaultdict(list)
        #print(i)
        for k,v in i:
            dictt[k].append(v)
            #print(dictt[k]) εάν το βάλω εδώ και το αρχικοποιήσω και στην αρχή τοτε εκτυπώνει ένα λεξικο με genotypes as keys και μία λίστα απ όλα τα dreqs του συγκεκριμένου genotype για όλες τις γενιές. 
        dictionaries_genot_freq_list.append(dict(dictt))
    #print(dictionaries_genot_freq_list)
    dictionaries_genot_freq_list_new=[] #list with three dictionaries, one for each generation,with genotypes as keys an float frequencies as values
    for i in dictionaries_genot_freq_list:
        dict1={}
        for k,v in i.items():
            for j in v:
                dict1[k]=float(j)
        dictionaries_genot_freq_list_new.append(dict1)
        #print(dict1)
    #print(dictionaries_genot_freq_list_new) 
    gen_dict_new=dict(zip(gen_list,dictionaries_genot_freq_list))# full dictionary with generations as keys and for each generation key a dictionary as values with genotypes:frequencies with not duplicate keys as genotypes.
    #print(gen_dict_new)
    gen_dict_final=dict(zip(gen_list,dictionaries_genot_freq_list_new)) #full dictionary with generations as keys and dictionaries of genotypes and float frequencies as values.
    print(gen_dict_final)
    gen_dict=dict(zip(gen_list,(genot_freq_new))) # dictionary generation i: list of genotypes:frequencies
    #print(gen_dict)
    half_freq_thres=(max(float_freq))/2
    enum_genot_freq_dict=dict(zip(genot_enum,float_freq)) # dictionary enumerated genotype:float frequency
    print(enum_genot_freq_dict)
    genot_dict=dict(zip(genot_enum,genot)) # dictionary enumerated genotype:str(genotype)
    #print(genot_dict)
    #matter_genots=[str(genot_dict[k])+":"+str(v) for k,v in enum_genot_freq_dict.items() if v>half_freq_thres] # genotypes that their frequency is higher than the threshold
    matter_genots=[genot_dict[k] for k,v in enum_genot_freq_dict.items() if v>half_freq_thres]
    matter_genots_unique=list(set(matter_genots))
    #print(matter_genots)
    print(matter_genots_unique)
    #for i in matter_genots:
        #print(i)
     #   for k,v in gen_dict.items():
      #      if i in v:
       #         print(k,i.split(":")[1])
       

    #first100pairs = {k: enum_genot_freq_dict[k] for k in list(enum_genot_freq_dict)[:100]}
   
    plt.bar(enum_genot_freq_dict.keys(),enum_genot_freq_dict.values(),width=100,color='g')
    #plt.bar(range(0,len(enum_genot_freq_dict),500),enum_genot_freq_dict.values(),width=10,color='g')
    #plt.bar(first100pairs.keys(),first2pairs.values(),width=10,color='g')
    plt.title("All generations Genotypes vs Frequencies plot")
    plt.xlabel("Enumerated genotypes")
    plt.ylabel("Frequencies")
    plt.show()
    
    df=pd.DataFrame(gen_dict_final)
    df=df.transpose()
    print(df)
    
    """
   # df.plot(figsize=(10,6), xticks=range(0, 3)).legend(title='', bbox_to_anchor=(1, 1))
    df.plot(figsize=(10,6), xticks=range(len(gen_list))).legend(title='', bbox_to_anchor=(1, 1))
    plt.title("All genotypes evolution")
    plt.xlabel("Generations")
    plt.ylabel("Frequencies")
    plt.show()
    
    """
    
    
    df_new=df.loc[:,('{}'.format(i) for i in matter_genots_unique)]
    #print(df_new)
    #df_new.plot(figsize=(10,6), xticks=range(0, 3)).legend(title='', bbox_to_anchor=(1, 1))
    df_new.plot(figsize=(6000,6), xticks=range(0,len(gen_list),1000)).legend(title='', bbox_to_anchor=(1, 1))
    plt.title("High frequency genotypes evolution")
    plt.xlabel("Generations")
    plt.ylabel("Frequencies")
    plt.show()
    