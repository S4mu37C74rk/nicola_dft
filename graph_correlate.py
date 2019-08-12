#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 08:49:44 2019

@author: sc2195
"""
import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy import stats

data = []
with open('all_data.csv') as f:
    reader = csv.reader(f)
    data = [np.asarray(row) for row in reader]
    data = np.asarray(data)
    f.close()
    
methods = data[0]
data = data.transpose()

x = 1
y = 2

x_arr = []
for num in data[x][1:]:
    try:
        x_arr.append(float(num))
    except:
        x_arr.append(None)
        
y_arr = []
for num in data[y][1:]:
    try:
        y_arr.append(float(num))
    except:
        y_arr.append(None)

noneindex = [n for n in range(len(x_arr)) if not isinstance(x_arr[n], float)]
for n in range(len(y_arr)):
    if not isinstance(y_arr[n], float) and not n in noneindex:
        noneindex.append(n)
        
x_arr = [x_arr[n] for n in range(len(x_arr)) if not n in noneindex]
y_arr = [y_arr[n] for n in range(len(y_arr)) if not n in noneindex]

for n in range(1, len(methods)):
    print('{:}: key {:}  '.format(methods[n], n), end='')
print('\n')

#calculate best fit
slope, intercept, r_value, p_value, std_err = stats.linregress(x_arr, y_arr)
predict_y = [intercept + slope * i for i in x_arr]

fig = plt.figure(figsize=(16,10))
plt.scatter(x_arr, y_arr)
plt.plot(x_arr, predict_y, color='black')

#format plot
plt.title('Correlation of MEP min between {:} and {:} methods'.format(data[x][0], data[y][0]))
plt.legend(['Best fit, R^2 = {:.3f}'.format(r_value)])
plt.grid(which='both')
plt.xlabel("{:} - MEP min".format(data[x][0]))
plt.ylabel("{:} - MEP min".format(data[y][0]))

#save fig and display
plt.savefig('graphs/{:}vs{:}.png'.format(data[x][0], data[y][0]))
plt.show()