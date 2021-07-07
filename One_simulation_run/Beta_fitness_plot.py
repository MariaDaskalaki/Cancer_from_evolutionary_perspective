#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 21:19:37 2020

@author: maria
"""

######################## Code for plotting the beta fitness site scenario ####################

import pandas as pd
import matplotlib.pyplot as plt

column_names=['Positions', 'Fitness_distribution']

data_constant=pd.read_csv("/home/maria/Diplomatiki_parousiasi/Beta_distribution.txt", sep='\t', lineterminator='\n', names=column_names, header=None)
data_constant.plot(x='Positions', y='Fitness_distribution')
plt.xlabel("Cells_positions")
plt.ylabel("Fitness_values")
plt.title("Beta_fitness_distribution")
plt.show()
