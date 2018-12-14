import sys as sys
import numpy as np

from polypy import Read as rd
from polypy import Density as Dens
from polypy import Write as wr
from polypy import Generic as ge


from scipy import stats
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev

def system_volume(data, timestep, output=None):  
    '''
    system volume - Calculate the volume at each timestep and return a volume as a function of time plot
    Parameters
    ----------
    data      : dictionary containing the atomic trajectories, lattice vectors, timesteps and number of atoms. 
    timestep  : Timestep of MD simulation  :    Float 
    output    : Output file name           :    String         :   Default: Volume.png
            
    Returns
    -------
    volume    : Volume at each timestep   :    1D numpy array
    time      : Time                      :    1D numpy array
    '''
    if output is None:
        filename = "Volume.png"
    
    volume = np.array([])
    time = np.array([])

    for i in range(0, data['timesteps']):
        volume = np.append(volume, (np.prod(data['lv'][i])))
        time = np.append(time, (i * timestep))
        
    wr.line_plot(time, volume, "Timestep", "System Volume (" r'$\AA$' ")", filename)

    return volume, time

def conductivity(plane, volume, diff, temperature):
    '''
    conductivity - Calculate the ionic conductivity 
    Parameters
    ----------
    plane          : Total number of charge carriers       : Integer
    volume         : lattice vectors                       : Float
    diff           : diffusion coefficient                 : Float
    temperature    : Temperature                           : Integer
              
    Returns
    -------
    conductivity   : Conductivity                          : Float     
    '''
    volume = volume * (10 ** -30)
    diff = diff * (10 ** -9)
    conc = plane / volume
    EV = ev ** 2
    constants = kb * temperature
    conductivity = ((diff * conc) * EV) / constants
    
    return conductivity

def three_d_diffusion_coefficient(x):
    '''
    Calculate the diffusion coefficient from the slope of MSD vs Time
    Parameters
    ----------
    x  : Gradient of 1D diffusion   : Float
    Return
    ------
    Overal Diffusion coefficient
    '''

    return ((np.average(x)) / 6) * 10

def one_d_diffusion_coefficient(x):
    '''
    Calculate the diffusion coefficient from the slope of MSD vs Time
    Parameters
    ----------
    x  : Gradient of 1D diffusion   : Float
    Return
    ------
    Overal Diffusion coefficient
    '''

    return ((np.average(x)) / 2) * 10
    
def linear_regression(x, y):
    '''
    msd_stats - Linear Regression 
    Parameters
    ----------
    x : X coordinates
    y : Y coordinates
    Return
    ------
    slope  : Overal gradient   : Float
    '''
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
    return slope 