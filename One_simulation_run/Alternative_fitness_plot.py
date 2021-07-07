#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 15:12:40 2021

@author: maria
"""

######################## Code for plotting the sin fitness site scenario #################

import pandas as pd
import matplotlib.pyplot as plt

column_names=['Cells_Positions', 'Fitness_distribution']
data=pd.read_csv("/home/maria/Diplomatiki_parousiasi/Alternative_fitness_site_distribution.txt", sep='\t', lineterminator='\n',names=column_names, header=None)
#data['Cells_Positions']=pd.to_numeric(data['Cells_Positions'])
#print(data)

data.plot(x='Cells_Positions', y='Fitness_distribution')
plt.xlabel('Cells_positions')
plt.ylabel("Fitness_values")
#plt.legend('Fitness_distribution', loc='upper right')
plt.title("Multiple_vertex_fitness_distribution")
plt.show()
