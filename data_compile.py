#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 16:45:31 2019

@author: sc2195
"""
import csv
from os import listdir
from os.path import isfile, join

KCALMOL2HARTREE = 627.509

with open('DFT_complete.csv') as f:
    reader = csv.reader(f)
    keys = [row[0] for row in reader]
    f.close()

onlyfiles = [f for f in listdir() if isfile(join(f))]
files = [i for i in onlyfiles if 'complete.csv' in i]

rows = []
rows.append([file[:file.index('_com')] for file in sorted(files)])
rows[0].insert(0, '')
rows[0].append('Experiment')

for key in keys:
    row = [key]
    for file in sorted(files):
        data = []
        with open(file) as f:
            reader = csv.reader(f)
            data = [row[2] for row in reader if row[0] == key]
            if 'ATX' in file:
                try: 
                    datum = '{:.9f}'.format(abs(float(data[0]))/KCALMOL2HARTREE)
                    row.append(datum)
                except:
                    row.append(None)
            else:
                try: 
                    datum = '{:.9f}'.format(abs(float(data[0])))
                    row.append(datum)
                except:
                    row.append(None)
        f.close()

    with open('nicola/_{:}.metadata'.format(key)) as meta:
        lineList = meta.readlines()
        betavals = [line for line in lineList if 'ExpBeta' in line]
        beta = float(betavals[0][15:].strip())
        row.append(beta)
        meta.close()

    rows.append(row)    

with open('all_data.csv', 'w') as db:
    writer = csv.writer(db)
    writer.writerows(rows)
