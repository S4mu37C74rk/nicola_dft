#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 16:01:37 2019

@author: sc2195
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.text import Annotation
import numpy as np
import csv

def clean_data(data, x, y):
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
    label_list = [data[0][n] for n in range(1, len(data[0])) if not n in noneindex]
    return x_arr, y_arr, label_list

# define two colors, just to enrich the example
labels_color_map = {0: '#20b2aa', 1: '#ff7373'}

# set the examples count
no_examples = 50

data = []
with open('all_data.csv') as f:
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
fig = plt.figure(figsize=(20, 16))
ax = plt.subplot()

# extract the scatterplot drawing in a separate function so we ca re-use the code
def draw_scatterplot():
    ax.scatter(
        axis_values_x,
        axis_values_y,
        picker=True
    )


# draw the initial scatterplot
draw_scatterplot()


# create and add an annotation object (a text label)
def annotate(axis, text, x, y):
    text_annotation = Annotation(text, xy=(x, y), xycoords='data')
    axis.add_artist(text_annotation)


# define the behaviour -> what happens when you pick a dot on the scatterplot by clicking close to it
def onpick(event):
    # step 1: take the index of the dot which was picked
    ind = event.ind

    # step 2: save the actual coordinates of the click, so we can position the text label properly
    label_pos_x = event.mouseevent.xdata
    label_pos_y = event.mouseevent.ydata

    # just in case two dots are very close, this offset will help the labels not appear one on top of each other
    offset = 0

    # if the dots are to close one to another, a list of dots clicked is returned by the matplotlib library
    for i in ind:
        # step 3: take the label for the corresponding instance of the data
        label = generated_labels[i]

        # step 5: create and add the text annotation to the scatterplot
        annotate(
            ax,
            label,
            label_pos_x + offset,
            label_pos_y + offset
        )

        # step 6: force re-draw
        ax.figure.canvas.draw_idle()

        # alter the offset just in case there are more than one dots affected by the click
        offset += 0.01


# connect the click handler function to the scatterplot
fig.canvas.mpl_connect('pick_event', onpick)

# create the "clear all" button, and place it somewhere on the screen
ax_clear_all = plt.axes([0.0, 0.0, 0.1, 0.05])
button_clear_all = Button(ax_clear_all, 'Clear all')


# define the "clear all" behaviour
def onclick(event):
    # step 1: we clear all artist object of the scatter plot
    ax.cla()

    # step 2: we re-populate the scatterplot only with the dots not the labels
    draw_scatterplot()

    # step 3: we force re-draw
    ax.figure.canvas.draw_idle()

# link the event handler function to the click event on the button
button_clear_all.on_clicked(onclick)

ax.title.set_text('Correlation of MEP min between {:} and {:} methods'.format(data[x][0], data[y][0]))
ax.grid(which='both')
ax.set_xlabel("{:} - absolute MEP min".format(data[x][0]))
ax.set_ylabel("{:} - absolute MEP min".format(data[y][0]))

# initial drawing of the scatterplot
plt.plot()
plt.show()
