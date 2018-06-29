import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt

import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta

from scipy import stats
from scipy.constants import codata

def pbc(r_new, r_old, Vec, c):
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
        
 #       if Cross == False and c == 1:
    #        print(r_new)
  #      elif Cross == True and c == 1:
 #           print("r_old - ", r_old, "r_new pre correction - ", r, "r_new - ", r_new, "correction - ", val, "coord - ", c, "vector -", Vec)
   
    else:
 #       if c == 1:
 #           print("using multi-box correction")
        r = r_new
        if (r_new - r_old) > Vec * 0.5:
            r_new = r_new - (Vec * val)                    
            Cross = True
        
        elif -(r_new - r_old) > Vec * 0.5:
            r_new = r_new + (Vec * val)  
            Cross = True
      #  if c == 1:
          #  print("r_old - ", r_old, "r_new pre correction - ", r, "r_new - ", r_new, "correction - ", val, "coord - ", c, "vector -", Vec)
    
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