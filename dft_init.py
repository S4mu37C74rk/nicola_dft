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
    data = []
    try:
        with open('results_{:}.csv'.format(basis), 'r') as db:
            reader = csv.reader(db)

            for row in reader:
                data.append(row)
            db.close()
    except:
        pass
    data.append(results)

    with open('results_{:}.csv'.format(basis), 'w') as db:
        writer = csv.writer(db)
        writer.writerows(data)
        db.close()
    pass

#load files and extract atom coordinates and the inchikey
onlyfiles = [f for f in listdir('nicola') if isfile(join('nicola', f))]
files = sorted([i for i in onlyfiles if '.pdb' in i])

search = input("Enter Inchi of required compound (enter to pass): ")
for file in files:
    if search == '':
        break
    elif search in file:
        LOGGER.info("Requested file at index {:}".format(files.index(file)))

basis = 'def2-TZVPP'
start = int(input("Enter index of first file: "))
end = int(input("Enter index of final file: "))

for file in files[start:end]:
    LOGGER.info("Molecule started {:}".format(time.ctime()))
    
    
    cycle_res = []
    for cycle in range(1):
        result = nicola_dft_reader.calc(file, basis, cycle)
        cycle_res.append(result)
    
    potmax = max([result[1] for result in cycle_res])
    potmin = min([result[2] for result in cycle_res])
    name = cycle_res[0][0]
    SaveResults([name, potmax, potmin], basis)
        
