#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 13:55:25 2019

@author: sc2195
"""
import matplotlib.pyplot as plt
import csv

types = {
'MOMOL_MNDO': 'b',
'MOMOL_MNDOD': 'r',
'MOMOL_AM1': 'g',
'ATX': 'y',
'MPC_AM1': 'c',
'MPC_PM6': 'm',
'MPC_MNDO': 'silver',
'MPC_PM7': 'lightgreen',
'MOMOL_PM7': 'navy',
'MOMOL_PM3': 'crimson',
'MPC_MNDOD': 'darkorange',
'MPC_PM3': 'lightpink',
'DFT': 'k',
'MOMOL_PM6': 'lime'
}

KCALMOL2HARTREE = 627.509
EV2HARTREE = 27

#form list of results files in directory
files = [
#'MOMOL_MNDO_complete.csv', 
#'MOMOL_MNDOD_complete.csv', 
#'MOMOL_AM1_complete.csv', 
#'ATX_complete.csv', 
#'MPC_AM1_complete.csv', 
#'MPC_PM6_complete.csv', 
#'MPC_MNDO_complete.csv', 
#'MPC_PM7_complete.csv', 
#'MOMOL_PM7_complete.csv', 
#'MOMOL_PM3_complete.csv', 
#'MPC_MNDOD_complete.csv', 
#'MPC_PM3_complete.csv', 
'DFT_complete.csv', 
#'MOMOL_PM6_complete.csv'
]

#initialise plot to map results to
fig = plt.figure(figsize=(16,10))
ax = fig.add_subplot(111)

for file in files:
    expt_beta = []
    MEP_min = []
    name = file.split('.')[0][:file.index('_com')]
    
    with open(file, 'r') as readFile:
        reader = csv.reader(readFile)
        for line in reader:
            if 'ATX' in file:
                expt_beta.append(float(line[4]))
                MEP_min.append(-1*float(line[2])/KCALMOL2HARTREE)
            elif 'DFT' in file:
                expt_beta.append(abs(float(line[4])))
                MEP_min.append(abs(float(line[2])))
            elif line[2] != '' and line[4] != '' and abs(float(line[2])) < 4:
                expt_beta.append(abs(float(line[4])))
                MEP_min.append(abs(float(line[2])))
                
    readFile.close()
    ax.scatter(expt_beta, MEP_min, s=10, c=types.get(name), label=name)

#format plot
plt.title('Expt_B against MEP min for all methods')
plt.grid(which='both')
plt.xlabel("Beta (exp)")
plt.ylabel("MEP min")
plt.legend(loc='upper left')

#save fig and display
plt.savefig('graphs/all.png')
plt.show()