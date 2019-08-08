#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 09:48:36 2019

@author: sc2195
"""
from pyscf import gto
from pyscf import dft
from pyscf import lib
from pyscf import df
from pyscf.tools import cubegen
import numpy as np
import csv
from os import listdir
from os.path import isfile, join

def dist(v1, v2=(0,0,0)):
    '''Returns the distance between two points in 3D space;
       will return distance from origin if no second vector is given.
    '''
    resv = np.asarray(v1)-np.asarray(v2)
    return np.sqrt(sum(resv**2))

#function takes isosurface array, mol format, and density matrix as arguments
def mep(surface, mol, dm):
    '''Calculate the MEP, by taking each point on the density isosurface.
       Returned as a NumPy array.
    '''
    coords = np.asarray([surface[i][1:] for i in range(len(surface))])
    
    #Nuclear potential at given points
    Vnuc = 0
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        Z = mol.atom_charge(i)
        rp = r - coords
        Vnuc += Z / np.einsum('xi,xi->x', rp, rp)**.5
    
    #Potential of electron density
    Vele = np.empty_like(Vnuc)
    for p0, p1 in lib.prange(0, Vele.size, 600):
        fakemol = gto.fakemol_for_charges(coords[p0:p1])
        ints = df.incore.aux_e2(mol, fakemol)
        Vele[p0:p1] = np.einsum('ijp,ij->p', ints, dm)
    
    MEP = Vnuc - Vele
    
    mep_arr = surface
    for i in range(len(mep_arr)):
        mep_arr[i][0] = MEP[i]
    
    return mep_arr

def getAtomCoordsfromPDB(file):
    '''Extract atom information from the PDB input file.
    '''
    with open(file) as f:
        lineList = f.readlines()
        data = [item for item in lineList if item[:6] == 'HETATM']
        atom_coords = ''
        for n in data:
            atom_coords += '     '+n.strip()[12:14]
            atom_coords += n.strip()[27:38]+','
            atom_coords += n.strip()[38:46]+','
            atom_coords += n.strip()[46:54]+';\n'
    f.close()    
        
    return(atom_coords)
    
def BuildandOptimise(atomcoordinates, basis):
    '''Build the atom representation, using the basis given. Dataset for
       Nicola is pre-optimised for high quality DFT calculations.
    '''
    molecule = gto.M(atom=atomcoordinates, basis=basis, unit='Angstrom')
    molecule.build()
    molecule.verbose = 0
    method = dft.RKS(molecule)
    method.grids.prune = dft.gen_grid.treutler_prune
    method.grids.atom_grid = {"H": (50, 194), "O": (50, 194),}
    method.xc = 'b3lypg'
    method.scf()
    
    return molecule, method

def ConvertCubetoArray(cube):
    '''Read density file .cube into Numpy array for processing.
    '''
    with open(cube) as f:
        
        lineList = f.readlines()
        natom = int(lineList[2][:6].strip())
        infoList = lineList[2:6+natom]
        dataList = lineList[6+natom:]
        
        dims = [item.split("   ")[1:] for item in infoList[1:4]]
        x, x_vec = int(dims[0][0]), np.asarray([float(i) for i in dims[0][1:]])
        y, y_vec = int(dims[1][0]), np.asarray([float(i) for i in dims[1][1:]])
        z, z_vec = int(dims[2][0]), np.asarray([float(i) for i in dims[2][1:]])
        dim_array = [[x, x_vec], [y, y_vec], [z, z_vec]]
        
        test_cube = []
        for line in dataList:
            row = line.split("  ")[1:]
            rowFlt = [float(item) for item in row]
            test_cube.extend(rowFlt)
        test_cube = np.split(np.asarray(test_cube),(x*y))
        test_cube = np.split(np.asarray(test_cube),x)
        test_cube = np.asarray(test_cube)
    f.close()
    
    return test_cube, dim_array

def CalcIsosurface(cube_data, dim_array, potential=0.002, tol=0.0005):
    '''Calculate the density isosurface from the .cube data.
    '''
    isosurface = np.empty([0,4])
    nvox = 0
    x, x_vec = dim_array[0][0], dim_array[0][1]
    y, y_vec = dim_array[1][0], dim_array[1][1]
    z, z_vec = dim_array[2][0], dim_array[2][1]
    
    #form isosurface Nx4 matrix with density at index 0 followed by x, y, and z coordinates
    for i in range(x):
        for j in range(y):
            for k in range(z):
                val = cube_data[i][j][k]
                if val < potential + tol and val > potential - tol:
                    vec = (i-x/2)*x_vec + (j-y/2)*y_vec + (k-z/2)*z_vec
                    vec = [np.insert(vec, 0, val, axis=0)]
                    isosurface = np.concatenate((isosurface, vec), axis=0)
                    nvox += 1
                elif val > potential - tol:
                    nvox += 1
    
    #calculate volume of molecule according to volume of voxel
    mol_volume = nvox*dist(x_vec)*dist(y_vec)*dist(z_vec)
    return isosurface, mol_volume

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

for file in files[:250]:
    atom_coords = getAtomCoordsfromPDB('nicola/'+file)
    key = file.split('.')[0][1:]
    print('Loading {:} ...'.format(file))
    
    basis = '6-311+g**'
    
    #optimise in chosen basis, calculate density in cube file, convert to array and then output chosen isosurface
    print('Building molecule in {:} basis'.format(basis.upper))
    molecule, method = BuildandOptimise(atom_coords, basis)
    cubegen.density(molecule, '{:}_den.cube'.format(key), method.make_rdm1())
    array, dim_array = ConvertCubetoArray('{:}_den.cube'.format(key))
    print('Generating isosurface ...  ', end='')
    isosurface, volume = CalcIsosurface(array, dim_array)
    print('Complete')
                    
    #call function that produces numpy array with potential at index 0, followed by x, y, and z coordinates
    print('Generating molecular electrostatic potential surface...  ', end='')
    mep_arr = mep(isosurface, molecule, method.make_rdm1())
    print('Complete')
    potentials = np.transpose(mep_arr)[0]
    potmax = max(potentials)
    potmin = min(potentials)
    results.append([key, potmax, potmin])
    print('Job complete:\n  InchiKey: {:}\n  Max. Pot: {:}\n  Min. Pot: {:}\n'.format(key, potmax, potmin))
    
SaveResults(results, basis)