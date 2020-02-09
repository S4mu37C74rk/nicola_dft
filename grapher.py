#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 16:01:37 2019

@author: sc2195
"""
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import operator
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures

plt.rcParams.update({'font.size': 20})

def clean_data(data, x, y):
    
    noneindex = []
    for n in range(len(data[x])):
        try:
            float(data[x][n])
            float(data[y][n])
        except:
            noneindex.append(n)
    
    x_arr, y_arr, label_list = [], [], []
    for n in range(1, len(data[x])):
        if not n in noneindex and not data[0][n] in label_list:
            x_arr.append(float(data[x][n]))
            y_arr.append(float(data[y][n]))
            label_list.append(data[0][n])
    
    return x_arr, y_arr, label_list

def regress(x, y):
    x = x[:, np.newaxis]
    x = np.insert(x, 0, [0], axis=0)
    y = y[:, np.newaxis]
    y = np.insert(y, 0, [0], axis=0)
    
    polynomial_features= PolynomialFeatures(degree=2, include_bias=False)
    x_poly = polynomial_features.fit_transform(x)
    
    model = LinearRegression(fit_intercept=False)
    model.fit(x_poly, y)
    y_poly_pred = model.predict(x_poly)
    
    r2 = r2_score(y,y_poly_pred)
    
    sort_axis = operator.itemgetter(0)
    sorted_zip = sorted(zip(x,y_poly_pred), key=sort_axis)
    x, y_poly_pred = zip(*sorted_zip)
    return x, y_poly_pred, r2

mode = input('Graph alpha or beta? (ENTER a OR b): ')
if mode == 'a':
    data_file = 'alpha_data.csv'
    graph_of = 'Alpha'
    search = 'max'
elif mode == 'b':
    data_file = 'all_data.csv'
    graph_of = 'Beta'
    search = 'min'
else:
    print('Invalid mode selected: {:}'.format(mode))

data = []
with open(data_file) as f:
    reader = csv.reader(f)
    data = [np.asarray(row) for row in reader]
    data = np.asarray(data)
    f.close()
    
methods = data[0]
data = data.transpose()

for n in range(1,len(methods)):
    print('{:} - key {:}  '.format(methods[n], n), end='')
print('\n')

#pick methods to graph
x = int(input('Select x axis data: '))
y = int(input('Select y axis data: '))
axis_values_x, axis_values_y, generated_labels = clean_data(data, x, y)


# draw a scatter-plot of the generated values
fig = plt.figure(figsize=(15, 10))
ax = plt.subplot()

# extract the scatterplot drawing in a separate function so we ca re-use the code
def draw_scatterplot():
    ax.scatter(
        axis_values_x,
        axis_values_y,
        picker=True
    )
    x, y_poly_pred, r2 = regress(np.asarray(axis_values_x), np.asarray(axis_values_y))
    ax.plot(x, y_poly_pred, color='black')
    return r2


# draw the initial scatterplot
r2 = draw_scatterplot()

ax.title.set_text('MEP {:} of {:} vs {:}'.format(search, data[y][0], data[x][0]))
ax.grid(which='both')
ax.legend(['Best fit, R^2 = {:.3f}'.format(r2)])
ax.set_ylabel("{:} - absolute MEP {:}".format(data[y][0], search))
ax.set_xlabel("{:} - absolute MEP {:}".format(data[x][0], search))

# initial drawing of the scatterplot
#plt.plot()
plt.savefig('report_graphs/{:}vs{:}_{:}.png'.format(methods[y], methods[x], graph_of))
#plt.show()
