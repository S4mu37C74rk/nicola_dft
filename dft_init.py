#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 16:17:52 2019

@author: sc2195
"""
import csv
import nicola_dft_reader
from os import listdir
from os.path import isfile, join

def SaveResults(results, basis):
    '''Save extracted results to a .csv file with the InchiKey. Suffixed with
       the basis used.
    '''
    with open('results_{:}.csv'.format(basis), 'w') as db:
          writer = csv.writer(db)
          writer.writerows(results)
    db.close()
    pass

results = []
#load files and extract atom coordinates and the inchikey
onlyfiles = [f for f in listdir('nicola') if isfile(join('nicola', f))]
files = [i for i in onlyfiles if '.pdb' in i]

basis = '6-311g**'
start = input("Enter index of first file: ")
end = input("Enter index of final file: ")

for file in files[start:end]:
    result = nicola_dft_reader.calc(file, basis)
    results.append(result)

SaveResults(results, basis)