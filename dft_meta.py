#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 01 16:28:26 2019

@author: sc2195
"""
import csv

    
lines = []

with open('nwchem_results.csv', 'r') as db:
     reader = csv.reader(db)
     for row in reader:
         alpha, beta, function = None, None, None
         try:
             with open('nicola/_{:}.metadata'.format(row[0]), 'r') as f:
                 lineList = f.readlines()
                 dataList = [i.split('     ') for i in lineList]
                 for d in dataList:
                     if d[0] == 'ExpBeta':
                         beta = float(d[1].strip())
                     elif d[0] == 'ExpAlpha':
                         alpha = float(d[1].strip())
                     elif 'OldClassType' in d[0]:
                         function = d[0].strip().split('  ')[2]
         except:
             pass
         
         row.extend([alpha, beta, function])
         lines.append(row)
     db.close()

with open('nwchem_results.csv', 'w') as writeFile:
     writer = csv.writer(writeFile)
     writer.writerows(lines)
     writeFile.close()