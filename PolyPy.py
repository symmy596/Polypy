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
Atom = "F"
Runs = 5
Bin = 5
Box = 0.1
UL = 15.0
LL = -5.0