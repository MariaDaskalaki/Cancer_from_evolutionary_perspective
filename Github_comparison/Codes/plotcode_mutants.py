#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 19:30:08 2020

@author: maria
"""

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import pandas as pd

with open ("Number of mutants.txt","r") as f:
    data=f.readlines()
    x=[line.split()[0]for line in data] #generation
    y=[line.split()[1]for line in data] #number of generations
    gen_list=[int(i) for i in y]
    mutants_number=[int(line.split()[2]) for line in data]
    final_dict=dict(zip(gen_list,mutants_number))
    #print(final_dict)
    #df=pd.DataFrame(final_dict)
    #print(df)
    

    random.seed(30)
    x_in_h=np.random.normal(size=10000-mutants_number[0]) # healthy
    y_in_h=np.random.normal(size=10000-mutants_number[0])
    z_in_h=np.random.normal(size=10000-mutants_number[0])

    random.seed(60)
    x_in_m=np.random.normal(size=mutants_number[0]) # mutants
    y_in_m=np.random.normal(size=mutants_number[0])
    z_in_m=np.random.normal(size=mutants_number[0])

    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    ax.scatter(x_in_h,y_in_h,z_in_h,c='b')
    ax.scatter(x_in_m,y_in_m,z_in_m,c='r')
    ax.legend(['healthy','mutants'])
    plt.title("Initial cells states")

    random.seed(30)
    x_m_h=np.random.normal(size=10000-mutants_number[int(len(gen_list)/2)]) # healthy
    y_m_h=np.random.normal(size=10000-mutants_number[int(len(gen_list)/2)])
    z_m_h=np.random.normal(size=10000-mutants_number[int(len(gen_list)/2)])

    random.seed(60)
    x_m_m=np.random.normal(size=mutants_number[int(len(gen_list)/2)]) # mutants
    y_m_m=np.random.normal(size=mutants_number[int(len(gen_list)/2)])
    z_m_m=np.random.normal(size=mutants_number[int(len(gen_list)/2)])

    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    ax.scatter(x_m_h,y_m_h,z_m_h,c='b')
    ax.scatter(x_m_m,y_m_m,z_m_m,c='r')
    ax.legend(['healthy','mutants'])
    plt.title("Middle cells states")

    random.seed(30)
    x_f_h=np.random.normal(size=10000-mutants_number[-1]) # healthy
    y_f_h=np.random.normal(size=10000-mutants_number[-1])
    z_f_h=np.random.normal(size=10000-mutants_number[-1])

    random.seed(60)
    x_f_m=np.random.normal(size=mutants_number[-1]) # mutants
    y_f_m=np.random.normal(size=mutants_number[-1])
    z_f_m=np.random.normal(size=mutants_number[-1])

    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    ax.scatter(x_f_h,y_f_h,z_f_h,c='b')
    ax.scatter(x_f_m,y_f_m,z_f_m,c='r')
    ax.legend(['healthy','mutants'])
    plt.title("Final cells states")
    plt.show()
    
    plt.plot(gen_list,mutants_number)
    #plt.xticks(np.arange(0,len(gen_list),100))
    plt.xlabel("Generations")
    plt.ylabel("Number of mutants")
    plt.title("Tumor size")
    plt.show()
    
