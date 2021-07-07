#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 00:03:27 2020

@author: maria
"""

#################### Code for plotting the linear fitness site scenario ####################3

import pandas as pd
import matplotlib.pyplot as plt

column_names=['Cells_Positions', 'Fitness_distribution']
data=pd.read_csv("/home/maria/Diplomatiki_parousiasi/Linear_fitness_site_distribution.txt", sep='\t', lineterminator='\n',names=column_names, header=None)
#data['Cells_Positions']=pd.to_numeric(data['Cells_Positions'])
#print(data)

data.plot(x='Cells_Positions', y='Fitness_distribution')
plt.xlabel('Cells_positions')
plt.ylabel("Fitness_values")
plt.title("Linear_fitness_distribution")
plt.show()

