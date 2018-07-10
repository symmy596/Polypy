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

sys.path.append('Path')
os.environ['QT_QPA_PLATFORM']='offscreen'

import time
start_time = time.time()


Atom = "Atom"
Bin = 0.1
box = 0.1
ul = 8.0
ll = 0.01
timestep = 0.25

# 1) Read Coordinates

natoms, timesteps, trajectories, lv = rd.read_history("HISTORY_F", Atom)

# 2) Total Number of Species within a given region

plane = ta.one_dimensional_density_sb(trajectories, ul=ul, ll=ll, direction="x")

# 3) One Dimensional Density Plot

ta.one_dimensional_density(trajectories, timesteps, lv, Bin=0.1)

# 4) Two Dimensional Density Plot

ta.two_dimensional_density(trajectories, timesteps, lv, box=box, direction='z', log=False)
