"""
MSD functions included with `polypy`. There are two MSD classes and one class to store the data generated from the MSD calculation.
The first class performs a standard MSD calculation for the entire dataset while the second class will perform an MSD calculation
within a specified region of the simulation cell.
"""

# Copyright (c) Adam R. Symington
# Distributed under the terms of the MIT License
# author: Adam R. Symington

import sys as sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from polypy import fig_params
from polypy import utils as ut
from polypy.read import Trajectory 
from scipy.constants import codata
from scipy import stats

class MSDContainer():
    """
    The :py:class:`polypy.msd.MSDContainer` class stores the output from the msd calculation.
    """
    def __init__(self):
        self.msd = np.array([])
        self.xymsd = np.array([])
        self.xzmsd = np.array([])
        self.yzmsd = np.array([])
        self.xmsd = np.array([])
        self.ymsd = np.array([])
        self.zmsd = np.array([])
        self.time = np.array([])
        self.sd = np.array([])
        self.sd_time = np.array([])

    def clean_data(self):
        """
        Post msd the data is a list of time vs msd for each run. This needs to be
        normalised to give one continuous series of points.
        """
        time, self.msd = self.smooth_msd_data(self.time, self.msd)
        self.xymsd = self.smooth_msd_data(self.time, self.xymsd)[1]
        self.xzmsd = self.smooth_msd_data(self.time, self.xzmsd)[1]
        self.yzmsd = self.smooth_msd_data(self.time, self.yzmsd)[1]
        self.xmsd = self.smooth_msd_data(self.time, self.xmsd)[1]
        self.ymsd = self.smooth_msd_data(self.time, self.ymsd)[1]
        self.zmsd = self.smooth_msd_data(self.time, self.zmsd)[1]
        self.time = time

    def smooth_msd_data(self, x, y):
        """
        Smooths the data from the smoothed msd function. The data consists
        of multiple msd runs but the data is unordered. This function will
        order the data and average all y values with equivalent x values.

        Args:
            x (:py:attr:`array_like`): Time data.
            y (:py:attr:`array_like`): MSD data.
        
        Returns:
            z (:py:attr:`array_like`): Time / MSD data.
        """
        xy = np.column_stack((x, y))
        z = pd.DataFrame(xy).groupby(0, as_index=False)[1].mean().values
        return z[:, 0], z[:, 1]

    def xyz_diffusion_coefficient(self, exclude_initial=0, exclude_final=1):
        """
        Calculates the three dimensional xyz diffusion coefficient.

        Returns:
            (:py:attr:`float`): xyx Diffusion coefficient.
        """
        gradient = stats.linregress(self.time[exclude_initial:-exclude_final], self.msd[exclude_initial:-exclude_final])[0]
        return (gradient / 6) * 10

    def xy_diffusion_coefficient(self, exclude_initial=0, exclude_final=1):
        """
        Calculates the two dimensional xy diffusion coefficient.

        Returns:
            (:py:attr:`float`): xy Diffusion coefficient.
        """
        gradient = stats.linregress(self.time[exclude_initial:-exclude_final], self.xymsd[exclude_initial:-exclude_final])[0]
        return (gradient / 4) * 10

    def xz_diffusion_coefficient(self, exclude_initial=0, exclude_final=1):
        """
        Calculates the two dimensional xz diffusion coefficient.

        Returns:
            (:py:attr:`float`): xz Diffusion coefficient.
        """
        gradient = stats.linregress(self.time[exclude_initial:-exclude_final], self.xzmsd[exclude_initial:-exclude_final])[0]
        return (gradient / 4) * 10

    def yz_diffusion_coefficient(self, exclude_initial=0, exclude_final=1):
        """
        Calculates the two dimensional yz diffusion coefficient.

        Returns:
            (:py:attr:`float`): yz Diffusion coefficient.
        """
        gradient = stats.linregress(self.time[exclude_initial:-exclude_final], self.yzmsd[exclude_initial:-exclude_final])[0]
        return (gradient / 4) * 10

    def x_diffusion_coefficient(self, exclude_initial=0, exclude_final=1):
        """
        Calculates the one dimensional x diffusion coefficient.

        Returns:
            (:py:attr:`float`): x Diffusion coefficient.
        """
        gradient = stats.linregress(self.time[exclude_initial:-exclude_final], self.xmsd[exclude_initial:-exclude_final])[0]
        return (gradient / 2) * 10
    
    def y_diffusion_coefficient(self, exclude_initial=0, exclude_final=1):
        """
        Calculates the one dimensional y diffusion coefficient.

        Returns:
            (:py:attr:`float`): y Diffusion coefficient.
        """
        gradient = stats.linregress(self.time[exclude_initial:-exclude_final], self.ymsd[exclude_initial:-exclude_final])[0]
        return (gradient / 2) * 10

    def z_diffusion_coefficient(self, exclude_initial=0, exclude_final=1):
        """
        Calculates the one dimensional z diffusion coefficient.

        Returns:
            (:py:attr:`float`): z Diffusion coefficient.
        """
        gradient = stats.linregress(self.time[exclude_initial:-exclude_final], self.zmsd[exclude_initial:-exclude_final])[0]
        return (gradient / 2) * 10


class MSD():
    """
    The :py:class:`polypy.msd.MSD` class calculates the mean squared displacements for a given atom.

    Args:
        data (:py:class:`polypy.read.Trajectory`): Object containing the information from the HISTORY or ARCHIVE files.
        sweeps (:py:attr:`int`, optional): How many times should the starting timestep be changed. Default is :py:attr:`1`. 
    """
    def __init__(self, data, sweeps=1):
        self.data = data
        self.sweeps = sweeps
        if self.data.timesteps == 1:
            raise ValueError("ERROR: - Only one timestep has been found")
        if len(np.unique(data.atom_name)) > 1:
            raise ValueError("ERROR: MSD can only handle one atom type. Exiting")
        if data.data_type == "DL_MONTE ARCHIVE":
            raise ValueError("DLMONTE simulations are not time resolved")
        self.distances = []
        self.msd_information = MSDContainer()

    def msd(self):
        """
        Calculates the mean squared displacement for the trajectory.

        Returns:
            (:py:class:`polypy.msd.MSDContainer`): Object containing the information for the MSD.
        """
        for i in range(1, self.sweeps+1):
            trajectories = np.split(self.data.fractional_trajectory, self.data.timesteps)
            distances, timestamp = self.calculate_distances(trajectories, i + 5)
            distances = np.asarray(distances)
            self.msd_information.time = np.append(self.msd_information.time, timestamp)
            self.squared_displacements(distances, i + 5)
        self.msd_information.clean_data()
        return self.msd_information

    def calculate_distances(self, trajectories, start):
        """
        Calculates the distances.

        Args:
            trajectories (:py:attr:`array_like`): Fractional coordinates.
            start (:py:attr:`float`): Timestep to start the calculation.

        Returns:
            distances (:py:attr:`array_like`): Distances.
            timestamp (:py:attr:`array_like`): Timesteps.
        """
        distances = []
        timestamp = []
        trajectories = np.asarray(trajectories)
        r0 = trajectories[start-1]
        rOd = trajectories[start-1]

        for j in range((start), self.data.timesteps):
            r1 = trajectories[j]
            fractional_distance = r1 - r0
            r1.tolist()
            rOd.tolist()

            for k in range(0, fractional_distance[:, 0].size):
                for i in range(0, 3):
                    cross, r_new = ut.pbc(r1[k, i], rOd[k, i])
                    if cross is True:
                        r1[k, i] = r_new
                        fractional_distance[k, i] = r_new - r0[k, i]
                cartesian_distance = np.matmul(self.data.lv[j], fractional_distance[k])
                distances.append(cartesian_distance)
            timestamp.append((j-start) * self.data.simulation_timestep)
            r1 = np.asarray(r1)
            rOd = np.asarray(rOd)
            rOd = r1
        return distances, timestamp

    def squared_displacements(self, distances, run):
        """
        Calculates the squared distances.

        Args:
            distances (:py:attr:`array_like`): Distances.
            run (:py:attr:`float`): Timestep to start the calculation.
        """
        squared_displacements = distances ** 2
        self.three_dimension_square_distance(squared_displacements, run)
        self.two_dimension_square_distance(squared_displacements, run)
        self.one_dimension_square_distance(squared_displacements, run)

    def three_dimension_square_distance(self, distances, run):
        """
        Calculate the MSD in three dimensions.

        Args:
            distances (:py:attr:`array_like`): Distances.
            run (:py:attr:`float`): Timestep to start the calculation.
        """
        summed_distances = np.sum(distances, axis=1)
        reshaped_array = np.reshape(summed_distances, (self.data.timesteps-run, self.data.total_atoms))
        msd = np.mean(reshaped_array, axis=1)
        self.msd_information.msd = np.append(self.msd_information.msd, msd)

    def two_dimension_square_distance(self, distances, run):
        """
        Calculate the MSD in two dimensions.

        Args:
            distances (:py:attr:`array_like`): Distances.
            run (:py:attr:`float`): Timestep to start the calculation.
        """
        xy_distances = np.sum(np.array([distances[:,0], distances[:,1]]), axis=0)
        xz_distances = np.sum(np.array([distances[:,0], distances[:,2]]), axis=0)
        yz_distances = np.sum(np.array([distances[:,1], distances[:,2]]), axis=0)
        xy_array = np.reshape(xy_distances, (self.data.timesteps-run, self.data.total_atoms))
        xz_array = np.reshape(xz_distances, (self.data.timesteps-run, self.data.total_atoms))
        yz_array = np.reshape(yz_distances, (self.data.timesteps-run, self.data.total_atoms))

        xy_msd = np.mean(xy_array, axis=1)
        xz_msd = np.mean(xz_array, axis=1)
        yz_msd = np.mean(yz_array, axis=1)

        self.msd_information.xymsd = np.append(self.msd_information.xymsd, xy_msd)
        self.msd_information.xzmsd = np.append(self.msd_information.xzmsd, xz_msd)
        self.msd_information.yzmsd = np.append(self.msd_information.yzmsd, yz_msd)

    def one_dimension_square_distance(self, distances, run):
        """
        Calculate the MSD in one dimension.

        Args:
            distances (:py:attr:`array_like`): Distances.
            run (:py:attr:`float`): Timestep to start the calculation.
        """
        x_array = np.reshape(distances[:,0], (self.data.timesteps-run, self.data.total_atoms))
        y_array = np.reshape(distances[:,1], (self.data.timesteps-run, self.data.total_atoms))
        z_array = np.reshape(distances[:,2], (self.data.timesteps-run, self.data.total_atoms))

        x_msd = np.mean(x_array, axis=1)
        y_msd = np.mean(y_array, axis=1)
        z_msd = np.mean(z_array, axis=1)

        self.msd_information.xmsd = np.append(self.msd_information.xmsd, x_msd)
        self.msd_information.ymsd = np.append(self.msd_information.ymsd, y_msd)
        self.msd_information.zmsd = np.append(self.msd_information.zmsd, z_msd)


class RegionalMSD():
    """
    The :py:class:`polypy.msd.RegionalMSD` class calculates the mean squared displacements for a given atom in a specific 
    region of a simulation cell. 

    Args:
        data (:py:class:`polypy.read.Trajectory`): Object containing the information from the HISTORY or ARCHIVE files.
        lower_boundary (:py:attr:`float`): Coordinate of the lower limit of the region of interest.
        upper_boundary (:py:attr:`int`, optional): Coordinate of the upper limit of the region of interest.
        dimension (:py:attr:`int`, optional): Direction perpedicular to the region of interest. Default is :py:attr:`'x'`.
        sweeps (:py:attr:`int`, optional): How many times should the starting timestep be changed. Default is :py:attr:`1`.
    """
    def __init__(self, data, lower_boundary, upper_boundary, dimension='x', sweeps=1, trajectory_length=100):
        self.data = data
        if self.data.timesteps == 1:
            raise ValueError("ERROR: - Only one timestep has been found")
        if len(np.unique(data.atom_name)) > 1:
            raise ValueError("ERROR: MSD can only handle one atom type. Exiting")
        if data.data_type == "DL_MONTE ARCHIVE":
            raise ValueError("DLMONTE simulations are not time resolved")
        self.lower_boundary = lower_boundary
        self.upper_boundary = upper_boundary
        self.trajectory_length = trajectory_length
        self.sweeps = sweeps
        if dimension == "x":
            self.dimension = 0
        elif dimension == "y":
            self.dimension = 1
        elif dimension == "z":
            self.dimension = 2

    def analyse_trajectory(self):
        """
        Analyse the trajectory object.

        Returns:
            msd_information (:py:class:`polypy.msd.MSDContainer`): MSDContainer object - MSD information.
        """
        xc = np.reshape(self.data.fractional_trajectory[:, self.dimension], ((self.data.timesteps),
                        self.data.total_atoms))

        trajectories = np.split(self.data.fractional_trajectory, self.data.timesteps)
        trajectories = np.asarray(trajectories)
        self.msd_information = MSDContainer()
        for i in range(0, self.data.total_atoms):
            self.check_trajectory(trajectories[:, i], xc[:, i])
        self.msd_information.clean_data()
        return self.msd_information

    def initialise_new_trajectory(self):
        """
        Create a new MSDContainer object, specific to slice of a trajectory.

        Returns:
            new_trajectory (:py:class:`polypy.msd.MSDContainer`): MSDContainer object.
        """
        new_trajectory = Trajectory(self.data.atom_list, "DLPOLY HISTORY")
        new_trajectory.record_number = self.data.record_number
        new_trajectory.time = self.data.time
        new_trajectory.total_atoms = 1
        new_trajectory.timesteps = self.data.timesteps
        new_trajectory.simulation_timestep = self.data.simulation_timestep
        return new_trajectory

    def update_msd_info(self, container):
        """
        Adds the information calculated for a single atom to the information 
        of the whole trajectory.

        Args:
            container (:py:class:`polypy.msd.MSDContainer`): MSDContainer object - single atom.
        """       
        self.msd_information.msd = np.append(self.msd_information.msd, container.msd)
        self.msd_information.xymsd = np.append(self.msd_information.xymsd, container.xymsd)
        self.msd_information.xzmsd = np.append(self.msd_information.xzmsd, container.xzmsd)
        self.msd_information.yzmsd = np.append(self.msd_information.yzmsd, container.yzmsd)
        self.msd_information.xmsd = np.append(self.msd_information.xmsd, container.xmsd)
        self.msd_information.ymsd = np.append(self.msd_information.ymsd, container.ymsd)
        self.msd_information.zmsd = np.append(self.msd_information.zmsd, container.zmsd)
        self.msd_information.time = np.append(self.msd_information.time, container.time)

    def check_trajectory(self, trajectory, xc):
        """
        Analyse the trajectory of an individual atom.

        Args:
            trajectory (:py:class:`polypy.read.Trajectory`): Trajectory object.
            xc (:py:attr:`array_like`): Coordinates perpendicular to the region of interest.
        """
        ib = False
        count = 0
        conductivity_count = 0
        new_trajectory = self.initialise_new_trajectory()
        for i in range(0, xc.size):
            if xc[i] > self.lower_boundary and xc[i] < self.upper_boundary:
                ib = True
                count = count + 1
                conductivity_count = conductivity_count + 1
                new_trajectory.fractional_trajectory.append(trajectory[i])
                new_trajectory.lv.append(self.data.lv[i])

            elif xc[i] < self.lower_boundary or xc[i] > self.upper_boundary:
                if count > self.trajectory_length and ib is True:
                    new_trajectory._clean_data()
                    atom_msd = MSD(new_trajectory, self.sweeps)
                    msd = atom_msd.msd()
                    self.update_msd_info(msd)
                    count = 0
                    new_trajectory = self.initialise_new_trajectory()
                else:

                    ib = False
                    new_trajectory = self.initialise_new_trajectory()
                    count = 0

        if count > self.trajectory_length and ib is True:
            new_trajectory._clean_data()
            atom_msd = MSD(new_trajectory, self.sweeps)
            msd = atom_msd.msd()
            self.update_msd_info(msd)
