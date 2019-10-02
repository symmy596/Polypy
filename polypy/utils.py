import numpy as np
from scipy import stats
from scipy.constants import codata
from scipy import integrate
import pandas as pd
from numpy import linalg as la

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev


def system_volume(data):
    '''Calculate the volume at each timestep and return a volume as a
    function of time plot.

    Parameters
    ----------


    Returns
    -------
    volume : array like
        Volume at each timestep
    time : array like
        Time
    '''
    volume = np.array([])
    step = np.array([])

    for i in range(data['timesteps']):
        volume = np.append(volume, (np.dot(data['lv'][i][0,:] , np.cross(data['lv'][i][1,:], data['lv'][i][2,:] ))))
        step = np.append(step, i)
    return volume, step


def conductivity(plane, volume, diff, temperature):
    '''Calculate the ionic conductivity

    Parameters
    ----------
    plane : int
        Total number of charge carriers
    volume : float
        System volume
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
    return slope, intercept, r_value, p_value, std_err


def pbc(rnew, rold, i):
    '''Periodic boundary conditions for an msd calculation

    Parameters
    ----------
    rnew : float
        Value of current atomic position
    rold : float
        Value of previous atomic position
    vec  : float
        Lattice vector at that timestep

    Return
    ------
    cross : bool
        Result of PBC check - True if atom crosses the boundary
    new : float
        New position
    '''
    shift = abs(rold - rnew)
  #  if i == 1:
  #      print(shift)
    shift = round(shift, 0)
  #  if i == 1:
  #      print(shift)
    shift = int(shift)
  #  if i == 1:
  #      print(shift)
    cross = False
  #  if i == 1:
  #      print("Old", rold, "New", rnew, "Org", r0, "timestep", j, "distance", distance)
    if shift < 2:
        if (rnew - rold) > 0.5:
            rnew = rnew - 1.0
            cross = True
        elif -(rnew - rold) > 0.5:
            rnew = rnew + 1.0
            cross = True
    else:
        if (rnew - rold) > 0.5:
            rnew = rnew - shift
            cross = True
        elif -(rnew - rold) > 0.5:
            rnew = rnew + shift 
            cross = True
   # if i == 1:
   #     print("New", rnew, "Old", rold, "Shift", shift, "New_Distance", (rnew-r0))

    return cross, rnew


def bin_choose(X, Y):
    '''Calculate the number of bins depending on a box size and a bin
    thickness

    Parameters
    ----------
    X : float
        box length
    Y : float
        bin thickness

    Returns
    -------
    Z  : float
        Number of bins
    '''
    Z = X / Y
    Z = round(Z, 0)
    Z = int(Z)
    Z = Z - 1
    return Z


def one_dimensional_charge_density(atoms_coords, atom_charges, bin_volume):
    """Calculates the charge density

    Parameters
    ----------
    atoms_coords : list
        List of numpy arrays containing number densities
        for a particlar species
    atom_charges : list
        List of charges corresponding to the atoms in atom_coords
    bin_volume : float
        Volume of the bins

    Returns
    -------
    charge_density : array like
        Charge density
    """
    number_density = np.column_stack((atoms_coords))
    charges = np.asarray(atom_charges)
    charge_density = np.sum(np.multiply(number_density, charges),
                            axis=1) / bin_volume
    return charge_density


def two_dimensional_charge_density(atoms_coords, atom_charges, bin_volume):
    """Calculates the charge density

    Parameters
    ----------
    atoms_coords : list
        List of numpy arrays containing number densities
        for a particlar species
    atom_charges : list
        List of charges corresponding to the atoms in atom_coords
    bin_volume : float
        Volume of the bins

    Returns
    -------
    charge_density : array like
        Charge density
    """
    number_density = np.dstack((atoms_coords))
    charges = np.asarray(atom_charges)
    charge_density = np.sum(np.multiply(number_density,
                                        charges), axis=2) / bin_volume
    return charge_density


def poisson_solver(dx, rho, timesteps):
    """Calculates the electric field and electrostatic potential from the
    charge density using the Poisson equation.

    Parameters
    ----------
    dx : array like
        Positions of each bin
    rho : array like
        Charge density
    timesteps :
        Total number of timesteps

    Returns
    -------
    dx : array like
        Positions of each bin
    e_field : array like
        Electric field
    potential : array like
        Electrostatic potential
    """
    magic_scaling = 14.3997584
    e_field = magic_scaling * integrate.cumtrapz(rho, dx, initial=0)
    e_field = e_field - np.mean(e_field)
    potential = -integrate.cumtrapz(e_field, dx, initial=0)
    potential = potential / timesteps
    return dx, e_field, potential


def smooth_msd_data(x, y):
    '''Smooths the data from the smoothed msd function. The data consists
    of multiple msd runs but the data is unordered. This function will
    order the data and average all y values with equivalent x values.

    Parameters
    ----------
    x : array like
        Time
    y : array like
        MSD

    Returns
    -------
    z : array like
        Smoothed Time and MSD values.
    '''
    xy = np.column_stack((x, y))
    z = pd.DataFrame(xy).groupby(0, as_index=False)[1].mean().values
    return z[:, 0], z[:, 1]


def calculate_rcplvs(lv):
    rcplvs = la.inv(np.transpose(lv))
    lengths = la.norm(lv, axis=1)
    return rcplvs, lengths


def cart_2_frac(coord, lengths, rcplvs):
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