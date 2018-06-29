import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt

import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta

from scipy import stats
from scipy.constants import codata

def pbc(r_new, r_old, vec):
    '''
    unfold_trajectory -  Imagine a scenario where an atom travels from box a to box a+1. In our simulation it has reapperared on 
    the other side of box a, but in reality it is in its current positon * the lenght of that lattice vector. Because in an MSD
    calculation we are interested in the distance between the current position and some reference position, we need to ensure that
    the new position remains in it "real" position. So we shift it with PBC. This algorythm calculates the difference between the 
    position at r_old and the new position r_new. If this distance is greater than the lenght of half the cell vector then the 
    particle has moved across a boundary and needs to be shifted back. So its new position is updated by shifting it back by the 
    lenght of the lattice vector * the box it should be in. 
    
    Parameters
    ----------
    first  : float
             Value of current atomic position
    second : float
             Value of previous atomic position
    third  : float
             Lattice vector at that timestep
    
    
    Return
    ------
    first  : str
             Result of PBC check - True if atom crosses the boundary
    second : float
             New position
    '''
    shift = abs((rold - rnew) / vec)
    shift = round(shift, 0)
    shift = int(shift)

    cross = False
    if shift < 2:

        if (rnew - rold) > vec * 0.5:
            rnew = rnew - vec                    
            cross = True
        
        elif -(rnew - rold) > vec * 0.5:
            rnew = rnew + vec  
            cross = True
         
    else:
        
        if (rnew - rold) > vec * 0.5:
            rnew = rnew - (vec * shift)                    
            cross = True
        
        elif -(rnew - rold) > vec * 0.5:
            rnew = rnew + (vec * shift)  
            cross = True
    
    return cross, rnew


def bin_choose(X, Y):
    '''
    BinChoose - Calculate the number of bins depending on a box size and a bin thickness
    
    Parameters
    ----------
    first  : float
             box length
    second : float
             bin thickness
             
    Return
    ------
    first : int
            Number of bins
            
    '''
    Z = X / Y
    Z = round(Z, 0)
    Z = int(Z)
    Z = Z - 1
    return Z