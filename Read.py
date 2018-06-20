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

     

def ReadHistory(File, Atom):
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
    
    '''
    if os.path.isfile(File):
        Coords = []
        NConfigs = 0
        History = open(File, 'r')
        Name = False
        Count = 0

        
        for line in History:
            if Name:
                Name = False
                Coords.append(line.split()) 
            if line[0] == Atom[0]:
                Name = True
                Count = Count + 1
            if line[0] =="t":
                NConfigs = NConfigs + 1
            
        Coords = np.asarray(Coords, dtype=float)

        NAtoms = Count / NConfigs
        NAtoms = int(NAtoms)

    else:
        print("File cannot be found")
        sys.exit(0)
    
    if NAtoms == 0:
        print("No Atoms of specified type exist within the selected HISTORY file")
        sys.exit(0)
        
    History.close() 
        
    return NAtoms, NConfigs, Coords



