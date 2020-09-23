"""
Analysis functions
"""

# Copyright (c) Adam R. Symington
# Distributed under the terms of the MIT License
# author: Adam R. Symington

import numpy as np
from scipy import integrate
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev


class OneDimensionalChargeDensity:
    """
    The :py:class:`polypy.analysis.OneDimensionalChargeDensity` 
    class converts one dimensional number densitie into
    the charge density, electric field and electrostatic potential.

    Args:
        histogram_positions (:py:attr:`array_like`): Histogram locations.
        atom_densities (:py:attr:`list`): List of histograms.
        atom_charges (:py:attr:`list`): List of atom charges.
        histogram_volume (:py:attr:`float`): Volume of the histograms.
        timesteps (:py:attr:`float`): Simulation timestep.
    """

    def __init__(self, histogram_positions, atom_densities, atom_charges, histogram_volume, timesteps):
        self.histogram_positions = histogram_positions
        self.atom_densities = atom_densities
        self.atom_charges = atom_charges
        self.histogram_volume = histogram_volume
        self.timesteps = timesteps
        self.scale = 14.3997584

    def calculate_charge_density(self):
        r"""
        Calculates the charge density in one dimension.

        Returns:
            charge_density (:py:attr:`array_like`): Charge density.
        """
        number_density = np.column_stack((self.atom_densities))
        charges = np.asarray(self.atom_charges)
        charge_density = np.sum(np.multiply(number_density, charges),
                            axis=1) / self.histogram_volume
        return self.histogram_positions, charge_density/self.timesteps

    def calculate_electric_field(self):
        r"""
        Calculates the electric field.

        Returns:
            e_field (:py:attr:`array_like`): Electric field.
        """
        rho = self.calculate_charge_density()[1]
        e_field = self.scale * integrate.cumtrapz(rho, self.histogram_positions, initial=0)
        e_field = e_field - np.mean(e_field)
        return self.histogram_positions, e_field

    def calculate_electrostatic_potential(self):
        r"""
        Calculates the electrostatic potential.

        Returns:
            potential (:py:attr:`array_like`): Electrostatic potential.
        """
        rho = self.calculate_charge_density()[1]
        e_field = self.scale * integrate.cumtrapz(rho, self.histogram_positions, initial=0)
        e_field = e_field - np.mean(e_field)
        potential = -integrate.cumtrapz(e_field, self.histogram_positions, initial=0)
        potential = potential / self.timesteps
        return self.histogram_positions, potential

def system_volume(data):
    """
    Calculate the volume at each timestep and return a volume as function of time.

    Args:
        data (:py:class:`polypy.read.Trajectory`): polypy Trajectory object.

    Returns:
        volume (:py:attr:`array_like`): Volume as a function of timestep.
        step (:py:attr:`array_like`): Timestep.
    """
    volume = []
    step = []

    for i in range(data.timesteps):
        volume.append((np.dot(data.lv[i][0,:] , np.cross(data.lv[i][1,:], data.lv[i][2,:] ))))
        step.append(i)
    return volume, step

def conductivity(charge_carriers, volume, diff, temperature, hr):
    """
    Calculate the ionic conductivity.

    Args:
        charge_carriers (:py:attr:`float`): Number of charge carriers.
        volume (:py:attr:`float`): Average cell volume.
        diff (:py:attr:`float`): Diffusion coefficient.
        temperature (:py:attr:`float`): Temperature.
        hr (:py:attr:`float`): Haven ratio.

    Returns:
        conductivity (:py:attr:`float`): Ionic conductivity.
    """
    volume = volume * (10 ** -24)
    diff = diff * (10 ** -8)
    conc = charge_carriers / volume
    EV = ev ** 2
    constants = kb * temperature
    conductivity = ((diff * conc) * EV) / constants
    return conductivity * hr

def two_dimensional_charge_density(atoms_coords, atom_charges, bin_volume, timesteps):
    """
    Calculates the charge density in two dimensions.

    Args:
        atoms_coords (:py:attr:`list`): List of atomic coordinates
        atom_charges (:py:attr:`list`): List of atomic charges
        bin_volume (:py:attr:`float`): Volume of histograms

    Returns:
        charge_density (:py:attr:`array_like`): Charge density.
    """
    number_density = np.dstack((atoms_coords))
    charges = np.asarray(atom_charges)
    charge_density = np.sum(np.multiply(number_density,
                                        charges), axis=2) / bin_volume
    return (charge_density / timesteps)
