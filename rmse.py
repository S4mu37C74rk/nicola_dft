#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 17:06:31 2019

@author: sc2195
"""
import csv
import numpy as np

def rmse(list1, list2):
    errors = []
    for i in range(len(list2)):
        if list2[i] != 'nan':
            try:
                errors.append((float(list2[i])-float(list1[i]))**2)
            except:
                pass
    mean = sum(errors)/len(errors)
    return np.sqrt(mean)

with open('all_data.csv') as data:
    reader = csv.reader(data)
    rows = [row for row in reader]
    data.close()
    
transrow = [[row[i] for row in rows] for i in range(len(rows[0]))]
dft = transrow[2]
compare = [item for item in transrow[1:] if item != transrow[2]]

for method in compare:
    result = rmse(dft[1:], method[1:])
    print('{:} vs {:} - RMSE = {:.3f}'.format(method[0], dft[0], result))
    