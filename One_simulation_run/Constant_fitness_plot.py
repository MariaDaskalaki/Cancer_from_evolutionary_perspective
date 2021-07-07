#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 00:15:58 2020

@author: maria
"""

################### Code for plotting the constant fitness site scenario ############################

import pandas as pd
import matplotlib.pyplot as plt

column_names=['Positions', 'Fitness_distribution']

data_constant=pd.read_csv("/home/maria/Diplomatiki_parousiasi/Constant_fitness_site_distribution.txt", sep='\t', lineterminator='\n', names=column_names, header=None)
data_constant.plot(x='Positions', y='Fitness_distribution')
plt.xlabel('Cells_positions')
plt.ylabel("Fitness_values")
plt.title("Constant_fitness_distribution")
plt.show()
