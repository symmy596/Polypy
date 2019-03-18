import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt

from polypy import Read as rd
from polypy import Density as Dens
from polypy import Utils as Ut
from polypy import Write as wr

from scipy import stats
from scipy.constants import codata

def pbc(rnew, rold, vec):
    '''Periodic boundary conditions for an msd calculation

    Parameters
    ----------
    rnew : float
        Value of current atomic position
    rold : float
        Value of previous atomic position
    vec  : float
        Lattice vector at that timestep
    
    Return
    ------
    cross : bool
        Result of PBC check - True if atom crosses the boundary
    new : float
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
    '''Calculate the number of bins depending on a box size and a bin
    thickness

    Parameters
    ----------
    X : float
        box length
    Y : float
        bin thickness
             
    Returns
    -------
    Z  : float
        Number of bins
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