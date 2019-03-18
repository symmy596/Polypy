import sys as sys
import numpy as np
from polypy import Read as rd
from polypy import Density as Dens
from polypy import Write as wr
from polypy import Generic as g
from scipy import stats
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev

def system_volume(data, timestep, output=None):  
    '''Calculate the volume at each timestep and return a volume as a
    function of time plot.
    
    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.    
    timestep : float
        Timestep of MD simulation
    output : str
        Output file name
            
    Returns
    -------
    volume : array like
        Volume at each timestep
    time : array like
        Time
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
    '''Calculate the ionic conductivity

    Parameters
    ----------
    plane : int
        Total number of charge carriers
    volume : float
        lattice vectors
    diff : float
        diffusion coefficient
    temperature : int
        Temperature
              
    Returns
    -------
    conductivity : float
        Conductivity    
    '''
    volume = volume * (10 ** -30)
    diff = diff * (10 ** -9)
    conc = plane / volume
    EV = ev ** 2
    constants = kb * temperature
    conductivity = ((diff * conc) * EV) / constants
    
    return conductivity

def three_d_diffusion_coefficient(x):
    '''Calculate the diffusion coefficient from the slope of MSD vs Time
    
    Parameters
    ----------
    x : float
        Gradient of 1D diffusion

    Returns
    ------
    float
        Overal Diffusion coefficient
    '''

    return ((np.average(x)) / 6) * 10

def one_d_diffusion_coefficient(x):
    '''Calculate the diffusion coefficient from the slope of MSD vs Time

    Parameters
    ----------
    x : float
        Gradient of 1D diffusion

    Returns
    -------
    float
        Overal Diffusion coefficient
    '''
    return ((np.average(x)) / 2) * 10
    
def linear_regression(x, y):
    '''Linear Regression 
    
    Parameters
    ----------
    x : array like
        X coordinates
    y : array like
        Y coordinates
    Return
    ------
    slope : float
        Overal gradient
    '''
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
    return slope 

def charge_density(atoms_coords, atom_charges):
  #  return 
    pass

def poisson_solver():
 #   return dx, e_field, potential
    pass