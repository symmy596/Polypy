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

sys.path.append('Path')
os.environ['QT_QPA_PLATFORM']='offscreen'

Atom = "BR"
Bin = 0.1
box = 0.1
ul = 8.0
ll = 0.01
timestep = 0.25

# 1) Read Coordinates

data = rd.read_history("HISTORY_F", Atom)

# 2) Volume Calculation

volume, t = ta.system_volume(lv, timesteps, timestep)

# 3) Density
dens = ta.Density(data)

# 4) Total Number of Species within a given region

plane = dens.one_dimensional_density_sb(ul=ul, ll=ll)
# 5) One Dimensional Density Plot

dens.one_dimensional_density(Bin=0.1)
# 6) Two Dimensional Density Plot

dens.two_dimensional_density(box=0.1)

# 6) Combined one and two dimensional density

dens.two_dimensional_density(box=0.1, log=True)
# 7) MSD Calculation with One Trajectory Sweep

ta.msd(data, timestep, conductivity=True, temperature=1500)
# 8) MSD Calculation with Multiple Trajectory Sweeps

ta.smooth_msd(data, timestep, runs=10, conductivity=True, temperature=1500)
# 9) MSD Calculation Within a Specific Region of the System

ta.plane_msd(data, timestep, runs=10, ul=ul, ll=ll, direction="x", conductivity=False)








