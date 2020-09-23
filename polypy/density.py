"""
Density functions included with `polypy`. The Density class will determine generate a three dimensional grid that stores the total number of times
that an atom spends within a xyz grid point during the simulation.
"""

# Copyright (c) Adam R. Symington
# Distributed under the terms of the MIT License
# author: Adam R. Symington

import numpy as np
from polypy import read as rd
from polypy import utils as ut


class Density:
    """
    The :py:class:`polypy.density.Density` class evaluates the positions of all atoms in the simulation.

    Args:
        data (:py:class:`polypy.read.Trajectory`): Object containing the information from the HISTORY or ARCHIVE files.
        histogram_size (:py:attr:`float`, optional): Specifies the spacing between histograms.
        atom (:py:attr:`str`, optional): Specifies the atom to calculate the density for. 
    """

    def __init__(self, data, histogram_size=0.1, atom=None):
        self.data = data
        self.histogram_size = histogram_size
        self.atom = atom
        if self.atom is None and len(self.data.atom_list) > 1:
            print("There are multiple different atoms in the trajectory, subsequent analysis may be incorrect")
        elif len(self.data.atom_list) > 1 and self.atom:
            self.data = self.data.get_atom(atom)
        self.lengths = self.data.cell_lengths
        self.x_lim = None
        self.y_lim = None
        self.z_lim = None
        self.x = None
        self.y = None
        self.z = None
        self.find_limits()
        self.coords_map = np.zeros((self.x_lim , self.y_lim , self.z_lim ))
        self.build_map()

    def find_limits(self):
        """
        Determine the upper and lower limits of the simulation cell in all three dimensions.
        """
        self.x_lim = (np.amax(self.lengths[:, 0]) / self.histogram_size).astype(int)
        self.y_lim = (np.amax(self.lengths[:, 1]) / self.histogram_size).astype(int)
        self.z_lim = (np.amax(self.lengths[:, 2]) / self.histogram_size).astype(int)
        self.x = (np.arange(0, self.x_lim)) * self.histogram_size
        self.y = (np.arange(0, self.y_lim)) * self.histogram_size
        self.z = (np.arange(0, self.z_lim)) * self.histogram_size

    def build_map(self):
        """
        Constructs three dimensional grid of histogram_size *  histogram_size * histogram_size.
        containing a count for how many times an atom passes through each histogram_size ** 3
        cube. 
        """
        positions = self.data.fractional_trajectory.tolist()

        for position in positions:
            
            self.update_map(position)

    def update_map(self, position):
        """
        Determines the specific location of a given atom and adds it to the corresponding 
        location in the three dimensional map of atomic positions.
        """
        xbox = (position[0] * self.x_lim).astype(int)
        ybox = (position[1] * self.y_lim).astype(int)
        zbox = (position[2] * self.z_lim).astype(int)
        self.coords_map[xbox, ybox, zbox] = self.coords_map[xbox, ybox, zbox] + 1

    def one_dimensional_density(self, direction="x"):
        """
        Calculate the particle density within one dimensional histograms of a structure.

        Args:
            direction (:py:attr:`str`): The dimension perpendicular to the histograms.

        Returns:
            x (:py:attr:`array_like`): Locations of histograms.
            y (:py:attr:`array_like`): Size of histograms.
            bin_volume (:py:attr:`float`): Volume of histograms.
        """
        if direction == "x":
            val = [1, 2]
            x = self.x
        elif direction == "y":
            val = [0, 2]
            x = self.y
        elif direction == "z":
            val = [0, 1]
            x = self.z
        tmp = np.sum(self.coords_map, axis=val[1])
        y = np.sum(tmp, axis=val[0])
        bin_volume = 0.1 * np.mean(self.data.cell_lengths[:,val[0]]) * np.mean(self.data.cell_lengths[:,val[1]])
        return x, y, bin_volume

    def two_dimensional_density(self, direction="x"):
        """
        Calculate the particle density within two dimensional pixels of a structure.

        Args:
            direction (:py:attr:`str`): The dimension normal to the pixels.

        Returns:
            x (:py:attr:`array_like`): Locations of one dimension of the pixels.
            y (:py:attr:`array_like`): Locations of one dimension of the pixels.
            z (:py:attr:`array_like`): Size of pixels.
            bin_volume (:py:attr:`float`): Volume of pixels.
        """
        if direction == "x":
            val = 0
            x = self.z
            y = self.y
        elif direction == "y":
            val = 1
            x = self.z
            y = self.x
        elif direction == "z":
            val = 2
            x = self.y
            y = self.x
        z = np.sum(self.coords_map, axis=val)
        box_volume = 0.1 * 0.1 * np.mean(self.data.cell_lengths[:,val])
        return x, y, z, box_volume
