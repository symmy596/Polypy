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

     

def read_history(file, atom):
    '''
    ReadHistory - Read a DL_POLY HISTORY file
    
    Parameters
    ----------
    file   : HISTORY filename                         : String
    atom   : Atom to be read                          : String
             
    Return
    ------
    data   : trajectories - Atomic trajectories       : Numpy array
           : lv           - Lattice Vectors           : Numpy array
           : timesteps    - Total number of timesteps : Integer
           : natoms       - Tital number of atoms     : Integer
    '''
    if os.path.isfile(file):
        trajectories = []
        timesteps = 0
        history = open(file, 'r')
        name = False
        count = 0
        tstep = False
        c = 0
        lv = []
        for line in history:
            if c == 3:
                c = 0
                tstep = False
            if c < 3 and tstep == True:
                lv.append(line.split())
                c = c + 1
            if name:
                name = False
                trajectories.append(line.split()) 
            if line[0] == atom[0]:
                name = True
                count = count + 1
            if line[0] =="t":
                timesteps = timesteps + 1
                tstep = True
        trajectories = np.asarray(trajectories, dtype=float)
        lv = np.asarray(lv, dtype=float)
        natoms = count / timesteps
        natoms = int(natoms)
        vec = np.array([])
        lv = np.split(lv, timesteps)

        for i in range(0, timesteps):
            vec = np.append(vec, (lv[i].sum(axis=0)))
        lv = np.reshape(vec, (timesteps, 3))

        data = {'trajectories':trajectories, 'lv':lv, 'timesteps':timesteps, 'natoms':natoms}
    
    else:
        print("File cannot be found")
        sys.exit(0)
    
    if natoms == 0:
        print("No Atoms of specified type exist within the selected HISTORY file")
        sys.exit(0)
        
    history.close() 
    return data
#


