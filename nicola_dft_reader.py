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
from pyscf import __config__
import os

setattr(__config__, 'cubegen_box_margin', 5.0)
from pyscf.tools import cubegen
import numpy as np
import logging

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

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

def getCoords(atom_info, cycle):
    '''Fetch coords, also performing a rotation according to the cycle
    '''
    x = float(atom_info[27:38].strip())
    y = float(atom_info[38:46].strip())
    z = float(atom_info[46:54].strip())
    vec = [x,y,z]
    
    return vec[(0+cycle)%3], vec[(1+cycle)%3], vec[(2+cycle)%3]

def getAtomCoordsfromPDB(file, cycle):
    '''Extract atom information from the PDB input file.
    '''
    with open(file) as f:
        lineList = f.readlines()
        data = [item for item in lineList if item[:6] == 'HETATM']
        
        atom_coords = ''
        for n in data:
            atom_info = n.strip()
            x, y, z = getCoords(atom_info, cycle)
            atom_coords += '     '+atom_info[12:14]+' '+str(x)+', '+str(y)+', '+str(z)+';\n'
            
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
    method.xc = 'b3lyp'
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
        
        dims = [[item[:7], item[7:19], item[19:31], item[31:]] for item in infoList[1:4]]
        dim_array = [[int(item[0].strip()), np.asarray([float(i.strip()) for i in item[1:]])] for item in dims]
        x, y = dim_array[0][0], dim_array[1][0]
        
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

def CalcIsosurface(cube_data, dim_array, potential=0.002, tol=0.00005):
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

def calc(file, basis, cycle):
    atom_coords = getAtomCoordsfromPDB('nicola/'+file, cycle)
    key = file.split('.')[0][1:]
    LOGGER.info('Loading {:} ...'.format(file))
    
    #optimise in chosen basis, calculate density in cube file, convert to array and then output chosen isosurface
    LOGGER.info('Building molecule in {:} basis'.format(basis.upper()))
    molecule, method = BuildandOptimise(atom_coords, basis)
    cubegen.density(molecule, '{:}_den.cube'.format(key), method.make_rdm1(), resolution=(1/6))
    array, dim_array = ConvertCubetoArray('{:}_den.cube'.format(key))
    os.remove('{:}_den.cube'.format(key))
    LOGGER.info('Generating isosurface ...  ')
    isosurface, volume = CalcIsosurface(array, dim_array)
    del array
    LOGGER.info('Complete')
                    
    #call function that produces numpy array with potential at index 0, followed by x, y, and z coordinates
    LOGGER.info('Generating molecular electrostatic potential surface...  ')
    LOGGER.debug(isosurface)
    mep_arr = mep(isosurface, molecule, method.make_rdm1())
    LOGGER.debug(mep_arr)
    del isosurface
    del molecule
    del method
    LOGGER.info('Complete')
    potentials = np.transpose(mep_arr)[0]
    potmax = max(potentials)
    potmin = min(potentials)
    result = [key, potmax, potmin]
    LOGGER.info('Job complete:\n  InchiKey: {:}\n  Max. Pot: {:}\n  Min. Pot: {:}\n'.format(key, potmax, potmin))
    
    return result
