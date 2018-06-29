import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt

import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta

from scipy import stats
from scipy.constants import codata

def pbc(r_new, r_old, Vec):
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
    check = False
    val = abs((r_old - r_new) / Vec)
    val = round(val, 0)
    val = int(val)

    Cross = False
    if val < 2:
        r = r_new

        if (r_new - r_old) > Vec * 0.5:
            r_new = r_new - Vec                    
            Cross = True
        
        elif -(r_new - r_old) > Vec * 0.5:
            r_new = r_new + Vec  
            Cross = True
         
    else:
        r = r_new
        if (r_new - r_old) > Vec * 0.5:
            r_new = r_new - (Vec * val)                    
            Cross = True
        
        elif -(r_new - r_old) > Vec * 0.5:
            r_new = r_new + (Vec * val)  
            Cross = True
    
    return Cross, r_new


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