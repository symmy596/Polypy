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
import Read_tes as rt

from scipy import stats
from scipy.constants import codata

os.environ['QT_QPA_PLATFORM']='offscreen'


import time
start_time = time.time()


Atom = "F"
Runs = 5
Bin = 0.1
Box = 0.1
UL = 8.0
LL = 0.0
timestep = 0.25

NAtoms, NConfigs, Coords, lv = rd.read_history("../HISTORY_F", Atom)

volume, time = ta.system_volume(lv, NConfigs, timestep)

#1 dimensional density plot

ta.one_dimensional_density(Coords, NAtoms, NConfigs, lv, Bin, "x")

#2 dimensional density plot

ta.two_dimensional_density(Coords, NAtoms, NConfigs, lv, Box, 'z')

#Single MSD run

ta.msd(Coords, NConfigs, NAtoms, timestep, lv)

#MSD within a specific region of the system

ta.plane_msd(Coords, NConfigs, NAtoms, UL, LL, "x", Runs, lv, timestep)

#Smoothed MSD

ta.smooth_msd(Coords, NConfigs, NAtoms, lv, timestep, Runs)

#MSD that plots diffusion coefficient of each atom against its average position 

ta.pmsd(Coords, lv, NConfigs, NAtoms, timestep, Bin, "x")


print("--- %s seconds ---" % (time.time() - start_time))





