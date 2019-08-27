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
import time
import logging

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

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
files = sorted([i for i in onlyfiles if '.pdb' in i])

search = input("Enter Inchi of required compound (enter to pass): ")
for file in files:
    if search == '':
        break
    elif search in file:
        LOGGER.info("Requested file at index {:}".format(files.index(file)))

basis = '6-311G**'
start = int(input("Enter index of first file: "))
end = int(input("Enter index of final file: "))

for file in files[start:end]:
    LOGGER.info("Molecule started {:}".format(time.ctime()))
    result = nicola_dft_reader.calc(file, basis, 0)
    result2 = nicola_dft_reader.calc(file, basis, 1)
    result3 = nicola_dft_reader.calc(file, basis, 2)
    results.append(result)
    results.append(result2)
    results.append(result3)

SaveResults(results, basis)
