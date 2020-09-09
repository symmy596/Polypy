"""
Util functions
"""

# Copyright (c) Adam R. Symington
# Distributed under the terms of the MIT License
# author: Adam R. Symington

import numpy as np
from scipy import stats
from scipy.constants import codata
from scipy import integrate
import pandas as pd
from numpy import linalg as la

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev

def pbc(rnew, rold):
    """
    Periodic boundary conditions for an msd calculation

    Args:
        rnew (:py:attr:`float`, optional): New atomic position
        rold (:py:attr:`float`, optional): Previous atomic position

    Returns:
        cross (:py:attr:`bool`, optional): Has the atom cross a PBC?
        rnew (:py:attr:`float`, optional): New position
    """
    shift = abs(rold - rnew)
    shift = round(shift, 0)
    shift = int(shift)
    cross = False
    if shift < 2:
        if rnew - rold > 0.5:
            rnew = rnew - 1.0
            cross = True
        elif -(rnew - rold) > 0.5:
            rnew = rnew + 1.0
            cross = True
    else:
        if rnew - rold > 0.5:
            rnew = rnew - shift
            cross = True
        elif -(rnew - rold) > 0.5:
            rnew = rnew + shift 
            cross = True
    return cross, rnew

def calculate_rcplvs(lv):
    """
    Convert cartesian lattice vectors to the fractional lattice vectors

    Args:
        lv (:py:attr:`array_like`, optional): Lattice vectors

    Returns:
        rcplvs (:py:attr:`array_like`, optional): Reciprcocal lattice vectors
        lengths (:py:attr:`array_like`, optional): Cell lengths
    """
    rcplvs = la.inv(np.transpose(lv))
    lengths = la.norm(lv, axis=1)
    return rcplvs, lengths

def cart_2_frac(coord, lengths, rcplvs):
    """
    Convert cartesian coordinates to the fractional coordinates

    Args:
        coord (:py:attr:`array_like`, optional): Cartesian coordinates
        lengths (:py:attr:`array_like`, optional): Cell lengths
        rcplvs (:py:attr:`array_like`, optional): Reciprcocal lattice vectors

    Returns:
        coords (:py:attr:`array_like`, optional): Reciprcocal coordinates
    """ 
    coords = []
    if coord.size > 3:
        for i in range(0, coord[:,0].size):
            frac = np.matmul( rcplvs, coord[i] )
            frac = np.mod( frac, 1 )
            coords.append(frac)
        coords = np.asarray(coords, dtype=float)
        coords = np.reshape(coords, (coord[:,0].size, 3))
    else:
        frac = np.matmul(rcplvs, coord)
        coords = np.mod(frac, 1)
    return coords