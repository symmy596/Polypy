import numpy as np
from polypy import write as wr
from scipy import stats
from scipy.constants import codata
from scipy import integrate

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
    return volume, time


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


def pbc(rnew, rold, vec):
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
    shift = abs((rold - rnew) / vec)
    shift = round(shift, 0)
    shift = int(shift)
    cross = False

    if shift < 2:
        if (rnew - rold) > vec * 0.5:
            rnew = rnew - vec
            cross = True
        elif -(rnew - rold) > vec * 0.5:
            rnew = rnew + vec
            cross = True
    else:
        if (rnew - rold) > vec * 0.5:
            rnew = rnew - (vec * shift)
            cross = True
        elif -(rnew - rold) > vec * 0.5:
            rnew = rnew + (vec * shift)
            cross = True

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


def get_integer(x, y):
    z = x / y
    z = round(z, 0)
    z = int(z)
    return z


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
    charge_density = np.sum(np.multiply(number_density, charges), axis=2) / bin_volume
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
