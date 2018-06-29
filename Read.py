import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math as mt

import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta

from scipy import stats
from scipy.constants import codata

     

def read_history(File, Atom):
    '''
    ReadHistory - Read a DL_POLY HISTORY file
    
    Parameters
    ----------
    first  : file
             HISTORY filename
    second : str
             Atom to be read
             
    Return
    ------
    first  : int
             Number of atoms
    second : int
             Number of timesteps
    third  : numpy array
             Atomic coordinates
    fourth : numpy array
             Numpy array containing the lattice vectors at each timestep
    
    '''
    if os.path.isfile(File):
        Coords = []
        NConfigs = 0
        History = open(File, 'r')
        Name = False
        Count = 0
        tstep = False
        c = 0
        lv = []
        for line in History:
            if c == 3:
                c = 0
                tstep = False
            if c < 3 and tstep == True:
                lv.append(line.split())
                c = c + 1
            if Name:
                Name = False
                Coords.append(line.split()) 
            if line[0] == Atom[0]:
                Name = True
                Count = Count + 1
            if line[0] =="t":
                NConfigs = NConfigs + 1
                tstep = True
        Coords = np.asarray(Coords, dtype=float)
        lv = np.asarray(lv, dtype=float)
        NAtoms = Count / NConfigs
        NAtoms = int(NAtoms)
        vec = np.array([])
        lv = np.split(lv, NConfigs)

        for i in range(0, NConfigs):
            vec = np.append(vec, (lv[i].sum(axis=0)))
        lv = np.reshape(vec, (NConfigs, 3))


    else:
        print("File cannot be found")
        sys.exit(0)
    
    if NAtoms == 0:
        print("No Atoms of specified type exist within the selected HISTORY file")
        sys.exit(0)
        
    History.close() 
        
    return NAtoms, NConfigs, Coords, lv



