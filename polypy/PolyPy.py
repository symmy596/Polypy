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


import time
start_time = time.time()


Atom = "BR"
Bin = 0.1
box = 0.1
ul = 8.0
ll = 0.01
timestep = 0.25

# 1) Read Coordinates

natoms, timesteps, trajectories, lv = rd.read_history("HISTORY_F", Atom)

# 2) Volume Calculation

volume, t = ta.system_volume(lv, timesteps, timestep)

# 3) Total Number of Species within a given region

plane = ta.one_dimensional_density_sb(trajectories, ul=ul, ll=ll, direction="x")

# 4) One Dimensional Density Plot

ta.one_dimensional_density(trajectories, timesteps, lv, Bin=0.1)

# 5) Two Dimensional Density Plot

ta.two_dimensional_density(trajectories, timesteps, lv, box=box, direction='z', log=False)

# 6) MSD Calculation with One Trajectory Sweep

ta.msd(trajectories, lv, timesteps, natoms, timestep, conductivity=True, temperature=1500)

# 7) MSD Calculation with Multiple Trajectory Sweeps

ta.smooth_msd(trajectories, lv, timesteps, natoms, timestep, runs=10, conductivity=True, temperature=1500)

# 8) MSD Calculation Within a Specific Region of the System

ta.plane_msd(trajectories, lv, timesteps, natoms, timestep, runs=10, ul=ul, ll=ll, direction="x", conductivity=False)


# 9) MSD Calculation - Returns Plot of Atomic Diffusion Coefficient vs Average Position 

ta.pmsd(trajectories, lv, timesteps, natoms, timestep)


print("--- %s seconds ---" % (time.time() - start_time))





