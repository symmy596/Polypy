import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt

import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta

from scipy import stats
from scipy.constants import codata

def pbc(rnew, rold, vec):
    '''
    pbc - Periodic boundary conditions for an msd calculation
    Parameters
    ----------
    rnew  : Value of current atomic position   : Float
    rold  : Value of previous atomic position  : Float
    vec   : Lattice vector at that timestep    : Float
    
    Return
    ------
    cross  : Result of PBC check - True if atom crosses the boundary   : Bool
    new    : New position                                              : Float
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
    X  : box length    : Float
    Y  : bin thickness : Float
             
    Return
    ------
    Z  : Number of bins : Float
    '''
    Z = X / Y
    Z = round(Z, 0)
    Z = int(Z)
    Z = Z - 1
    return Z

def get_integer(x, y):
    z = x / y
    z = round(z, 0)
    z = int(z)
    return z