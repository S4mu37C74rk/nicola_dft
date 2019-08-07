#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 09:25:09 2019

@author: sc2195
"""
from pyscf import gto
from pyscf import dft
from pyscf import lib
from pyscf import df
from pyscf.geomopt import berny_solver
from pyscf.tools import cubegen
import numpy as np

#returns the distance between two points in 3d space, takes 3-member lists as arg
def dist(v1, v2=(0,0,0)):
    resv = np.asarray(v1)-np.asarray(v2)
    return np.sqrt(sum(resv**2))

#function takes isosurface array, mol format, and density matrix as arguments
def mep(surface, mol, dm):
    
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

h2o = gto.M(atom='      O 0.00000000,  0.000000,  0.119748;\
                        H 0.00000000,  0.761561, -0.478993;\
                        H 0.00000000, -0.761561, -0.478993',\
                        basis='6-311+G**')

h2o.build()
h2o.verbose = 0

method = dft.RKS(h2o)
method.grids.prune = dft.gen_grid.treutler_prune
method.grids.atom_grid = {"H": (50, 194), "O": (50, 194),}
method.xc = 'b3lypg'
method.scf()

#optimise the geometry of the molecule
opt_mol = berny_solver.optimize(method, verbose=0)

#calculate electron density of the molecule and write to Gaussian cube file
cubegen.density(opt_mol, 'water_den.cube', method.make_rdm1())

#read density file into Numpy array for processing
with open('water_den.cube') as f:
    
    lineList = f.readlines()
    infoList = lineList[2:9]
    dataList = lineList[9:]
    
    dims = [item.split("   ")[1:] for item in infoList[1:4]]
    x, x_vec = int(dims[0][0]), np.asarray([float(i) for i in dims[0][1:]])
    y, y_vec = int(dims[1][0]), np.asarray([float(i) for i in dims[1][1:]])
    z, z_vec = int(dims[2][0]), np.asarray([float(i) for i in dims[2][1:]])
    
    test_cube = list()
    for line in dataList:
        row = line.split("  ")[1:]
        rowFlt = [float(item) for item in row]
        test_cube.extend(rowFlt)
    test_cube = np.split(np.asarray(test_cube),(x*y))
    test_cube = np.split(np.asarray(test_cube),x)
    test_cube = np.asarray(test_cube)

#define potential on isosurface and tolerance
potential = 0.002
tol = 0.0005
isosurface = np.empty([0,4])
nvox = 0

#form isosurface Nx4 matrix with density at index 0 followed by x, y, and z coordinates
for i in range(x):
    for j in range(y):
        for k in range(z):
            val = test_cube[i][j][k]
            if val < potential + tol and val > potential - tol:
                vec = (i-x/2)*x_vec + (j-y/2)*y_vec + (k-z/2)*z_vec
                vec = [np.insert(vec, 0, val, axis=0)]
                isosurface = np.concatenate((isosurface, vec), axis=0)
                nvox += 1
            elif val > potential - tol:
                nvox += 1

#calculate volume of molecule according to volume of voxel
vol = nvox*dist(x_vec)*dist(y_vec)*dist(z_vec)
                
#call function that produces numpy array with potential at index 0, followed by x, y, and z coordinates
mep_arr = mep(isosurface, opt_mol, method.make_rdm1())
potentials = np.transpose(mep_arr)[0]
print(min(potentials))
