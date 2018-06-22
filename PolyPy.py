import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math as mt


import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta
import Write as wr

from scipy import stats
from scipy.constants import codata

os.environ['QT_QPA_PLATFORM']='offscreen'


import time
start_time = time.time()

Vec = np.array([27.32, 27.32, 27.32])
Atom = "BR"
Runs = 5
Bin = 0.1
Box = 0.1
UL = 10.0
LL = 0.0

#1 dimensional density plot

ta.one_dimensional_density(Coords, NAtoms, NConfigs, Vec, Bin, "x")

#2 dimensional density plot

ta.two_dimensional_density(Coords, NAtoms, NConfigs, Vec, Box, 'z')

#Single MSD run

ta.msd(Coords, Vec, NConfigs, NAtoms)

#MSD within a specific region of the system
#- This needs work

ta.plane_msd(Coords, NConfigs, NAtoms, UL, LL, Vec)

#Smoothed MSD

ta.smooth_msd(Coords, Vec, Runs, NConfigs, NAtoms)

#MSD that plots diffusion coefficient of each atom against its average position 

ta.pmsd(Coords, Vec, NConfigs, NAtoms, Bin)


print("--- %s seconds ---" % (time.time() - start_time))
