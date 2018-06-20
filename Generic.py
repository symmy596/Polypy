import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt

import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta

from scipy import stats
from scipy.constants import codata

def PBC(r_new, r_old, Vec):
    '''
    PBC - Periodic boundary conditions
    
    Parameters
    ----------
    first  : float
             Value of current atomic position
    second : float
             Value of previous atomic position
    third  : float
             Lattice vector
    
    Return
    ------
    first  : str
             Result of PBC check - True if atom crosses the boundary
    second : float
             New position
    '''
    
    Cross = False
    if (r_new - r_old) > Vec * 0.5:
        r_new = r_new - Vec                    
        Cross = True
        
    elif -(r_new - r_old) > Vec * 0.5:
        r_new = r_new + Vec  
        Cross = True
    
    return Cross, r_new

def BinChoose(X, Y):
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