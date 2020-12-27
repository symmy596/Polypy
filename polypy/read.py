"""
Read functions of `polypy`. Herein contains classes to read DL_POLY HISTORY
/ CONFIG files and DL_MONTE ARCHIVE files. All of the data that is
extracted from these files is stored in a trajectory class that is compatible
with all three file types.
"""

# Copyright (c) Adam R. Symington
# Distributed under the terms of the MIT License
# author: Adam R. Symington

import os as os
import numpy as np
from polypy import utils as ut


class Trajectory():
    """
    The :py:class:`polypy.read.Trajectory` class evaluates the positions
    of all atoms in the simulation.

    Args:
        atom_list (:py:class:`list`): List of unique atom names in trajectory.
        datatype (:py:attr:`str`): Datatype of the original dataset
        e.g. DL_POLY HISTORY or CONFIG.
    """
    def __init__(self, atom_list, datatype):
        self.atom_list = atom_list
        self.data_type = datatype
        self.cartesian_trajectory = []
        self.fractional_trajectory = []
        self.reciprocal_lv = []
        self.atom_name = []
        self.lv = []
        self.cell_lengths = []
        self.atoms_in_history = 0
        self.timesteps = 0
        self.total_atoms = 0
        self.atoms_at_timestep = []
        self.record_number = []
        self.time = []
        self.simulation_timestep = None

    def _clean_data(self):
        """
        Converts the data from lists / floats to numpy arrays / integers.
        """
        self.cartesian_trajectory = np.asarray(self.cartesian_trajectory,
                                               dtype=float)
        self.fractional_trajectory = np.asarray(self.fractional_trajectory,
                                                dtype=float)
        self.atom_name = np.asarray(self.atom_name, dtype=str)
        self.lv = np.asarray(self.lv, dtype=float)
        self.reciprocal_lv = np.asarray(self.reciprocal_lv, dtype=float)
        self.cell_lengths = np.asarray(self.cell_lengths, dtype=float)
        self.timesteps = int(self.timesteps)
        self.total_atoms = int(self.total_atoms)
        if self.data_type == "DL_POLY HISTORY":
            self.simulation_timestep = (self.record_number[1] -
                                        self.record_number[0]) * self.time[0]
        
    def get_config(self, timestep):
        """
        Isolates a specific DL_POLY CONFIG from a HISTORY file.

        Args:
            timestep (:py:class:`int`): Timestep of desired CONFIG.

        Returns:
            config_trajectory (:py:class:`polypy.read.Trajectory`): 
            Trajectory object for desired CONFIG.
        """
        if self.data_type == "DL_POLY CONFIG":
            raise ValueError("Only one timestep was found")
        config_trajectory = Trajectory(self.atom_list, "DL_POLY CONFIG")
        config_trajectory.cartesian_trajectory = np.split(
                                                 self.cartesian_trajectory,
                                                 self.timesteps)[timestep]
        config_trajectory.fractional_trajectory = np.split(
                                                  self.fractional_trajectory,
                                                  self.timesteps)[timestep]
        config_trajectory.atom_name = np.split(self.atom_name,
                                               self.timesteps)[timestep]
        config_trajectory.reciprocal_lv = self.reciprocal_lv[timestep]
        config_trajectory.lv = self.lv[timestep]
        config_trajectory.cell_lengths = np.split(self.cell_lengths,
                                                  self.timesteps)[timestep]
        config_trajectory.atoms_in_history = self.total_atoms
        config_trajectory.record_number = self.record_number[timestep]
        config_trajectory.time = self.time[timestep]
        config_trajectory.simulation_timestep = self.simulation_timestep
        config_trajectory.timesteps = 1
        config_trajectory.total_atoms = self.total_atoms
        return config_trajectory

    def get_atom(self, atom):
        """
        Isolates the trajectory for a specific atom type.

        Args:
            atom (:py:class:`str`): Atom label.

        Returns:
            atom_trajectory (:py:class:`polypy.read.Trajectory`):
            Trajectory object for desired atom.
        """
        if np.unique(self.atom_list).size == 1:
            raise ValueError("Only one atom type was found")
        atom_trajectory = Trajectory([atom], self.data_type)
        for i in range(0, self.atom_name.size):
            if self.atom_name[i] == atom:
                atom_trajectory.cartesian_trajectory.append(
                                            self.cartesian_trajectory[i])
                atom_trajectory.fractional_trajectory.append(
                                            self.fractional_trajectory[i])
                atom_trajectory.atom_name.append(self.atom_name[i])
                atom_trajectory.atoms_in_history = atom_trajectory.atoms_in_history + 1
        atom_trajectory.reciprocal_lv = self.reciprocal_lv
        atom_trajectory.lv = self.lv
        atom_trajectory.cell_lengths = self.cell_lengths
        atom_trajectory.timesteps = self.timesteps
        atom_trajectory.total_atoms = atom_trajectory.atoms_in_history / self.timesteps
        atom_trajectory.time = self.time
        atom_trajectory.record_number = self.record_number
        atom_trajectory.timesteps = 1
        atom_trajectory._clean_data()
        return atom_trajectory

    def remove_initial_timesteps(self, timesteps_to_exclude):
        """
        Removes timesteps from the beggining of a simulation

        Args:
            timesteps_to_exclude (:py:class:`int`): Number of timesteps to exclude

        Returns:
            new_trajectory (:py:class:`polypy.read.Trajectory`): 
            Trajectory object.
        """
        rows_to_exclude = timesteps_to_exclude * self.total_atoms
        if self.data_type == "DL_POLY CONFIG":
            raise ValueError("Only one timestep was found")
        new_trajectory = Trajectory(self.atom_list, self.data_type)
        new_trajectory.cartesian_trajectory = self.cartesian_trajectory[rows_to_exclude:]
        new_trajectory.fractional_trajectory = self.fractional_trajectory[rows_to_exclude:]
        new_trajectory.atom_name = self.atom_name
        new_trajectory.reciprocal_lv = self.reciprocal_lv[timesteps_to_exclude:]
        new_trajectory.lv = self.lv[timesteps_to_exclude:]
        new_trajectory.cell_lengths = self.cell_lengths[timesteps_to_exclude:]
        new_trajectory.atoms_in_history = self.atoms_in_history - (self.total_atoms * timesteps_to_exclude)
        new_trajectory.record_number = self.record_number[timesteps_to_exclude:]
        new_trajectory.time = self.time[timesteps_to_exclude:]
        new_trajectory.simulation_timestep = self.simulation_timestep
        new_trajectory.timesteps = self.timesteps - timesteps_to_exclude
        new_trajectory.total_atoms = self.total_atoms
        return new_trajectory

    def remove_final_timesteps(self, timesteps_to_exclude):
        """
        Removes timesteps from the end of a simulation

        Args:
            timesteps_to_exclude (:py:class:`int`): Number of timesteps to exclude

        Returns:
            new_trajectory (:py:class:`polypy.read.Trajectory`): 
            Trajectory object.
        """
        rows_to_exclude = timesteps_to_exclude * self.total_atoms
        if self.data_type == "DL_POLY CONFIG":
            raise ValueError("Only one timestep was found")
        new_trajectory = Trajectory(self.atom_list, self.data_type)
        new_trajectory.cartesian_trajectory = self.cartesian_trajectory[:rows_to_exclude]
        new_trajectory.fractional_trajectory = self.fractional_trajectory[:rows_to_exclude]
        new_trajectory.atom_name = self.atom_name
        new_trajectory.reciprocal_lv = self.reciprocal_lv[:timesteps_to_exclude]
        new_trajectory.lv = self.lv[:timesteps_to_exclude]
        new_trajectory.cell_lengths = self.cell_lengths[:timesteps_to_exclude]
        new_trajectory.atoms_in_history = self.atoms_in_history - (self.total_atoms * timesteps_to_exclude)
        new_trajectory.record_number = self.record_number[:timesteps_to_exclude]
        new_trajectory.time = self.time[:timesteps_to_exclude]
        new_trajectory.simulation_timestep = self.simulation_timestep
        new_trajectory.timesteps = self.timesteps - timesteps_to_exclude
        new_trajectory.total_atoms = self.total_atoms
        return new_trajectory

class History():
    """
    The :py:class:`polypy.read.Trajectory` class evaluates the positions of all atoms in the simulation.

    Args:
        atom_list (:py:class:`list`): List of unique atom names in trajectory.
        datatype (:py:attr:`str`): Datatype of the original dataset e.g. DL_POLY HISTORY.
    """
    def __init__(self, file, atom_list):
        self.file = file
        if os.path.isfile(self.file) is False:
            raise ValueError("File does not exist")
        self.atom_list = atom_list
        self.data_type = "DL_POLY HISTORY"
        self.trajectory = Trajectory(self.atom_list,
                                     self.data_type)

        self.read_history()
        self.trajectory.total_atoms = self.trajectory.atoms_in_history / self.trajectory.timesteps
        self._check_data()
        self.trajectory._clean_data()

    def read_history(self):
        """
        Reads a DL_POLY HISTORY file line by line and updates a :py:class:`polypy.read.Trajectory` object.
        """
        c = 0
        name = False
        tstep = False
        current_lv = []
        history = open(self.file, 'r')
        for line in history:
            split_line = line.split()
            if c == 3:
                c = 0
                tstep = False
                current_lv = np.asarray(current_lv, dtype=float)
                rcplvs, length = ut.calculate_rcplvs(current_lv)
                self.trajectory.cell_lengths.append(length)
                self.trajectory.lv.append(current_lv)
                self.trajectory.reciprocal_lv.append(rcplvs)
                current_lv = []
            if c < 3 and tstep is True:
                if c == 4:
                    current_lv.insert(0, split_line)
                else:
                    current_lv.append(split_line)
                c = c + 1
            if name:
                name = False
                self.trajectory.cartesian_trajectory.append(split_line)
                frac = np.matmul(rcplvs, np.asarray(split_line, dtype=float))
                frac = np.mod(frac, 1)
                self.trajectory.fractional_trajectory.append(frac)
            if split_line[0] in self.atom_list:
                self.trajectory.atom_name.append(split_line[0])
                name = True
                self.trajectory.atoms_in_history = self.trajectory.atoms_in_history + 1
            if split_line[0] == "timestep":
                self.trajectory.timesteps = self.trajectory.timesteps + 1
                tstep = True
                self.trajectory.record_number.append(float(split_line[1]))
                self.trajectory.time.append(float(split_line[5]))
        history.close()

    def _check_data(self):
        """
        Error handling function.
        """
        if self.trajectory.total_atoms == int(self.trajectory.total_atoms) is False:
            raise ValueError("The total number of atoms is not constant across each timestep")
        if self.trajectory.total_atoms == 0:
            raise ValueError("No Atoms of specified type exist within the HISTORY file")


class Config:
    """
    The :py:class:`polypy.read.Trajectory` class evaluates the positions of all atoms in a CONFIG.

    Args:
        atom_list (:py:class:`list`): List of unique atom names in trajectory.
        datatype (:py:attr:`str`): Datatype of the original dataset e.g. DL_POLY CONFIG.
    """
    def __init__(self, file, atom_list):
        self.file = file
        if os.path.isfile(self.file) is False:
            raise ValueError("File does not exist")
        self.atom_list = atom_list
        self.data_type = "DL_POLY CONFIG"
        self.trajectory = Trajectory(self.atom_list,
                                     self.data_type)
        self.read_config()
        self.trajectory.timesteps = 1
        self._check_data()
        self.trajectory._clean_data()

    def read_config(self):
        """
        Read a DL_POLY HISTORY file line by line and updates a :py:class:`polypy.read.Trajectory` object.
        """
        config = open(self.file, 'r')
        name = False
        title = config.readline()
        stuff = config.readline()
        lv = []
        for i in range(0, 3):
            lline = config.readline()
            lv.append(lline.split())
        lv = np.asarray(lv, dtype=float)
        rcplvs, length = ut.calculate_rcplvs(lv)
        self.trajectory.cell_lengths.append(length)
        self.trajectory.reciprocal_lv.append(rcplvs)
        for line in config:
            split_line = line.split()
            if name:
                name = False
                self.trajectory.cartesian_trajectory.append(split_line)
                frac = np.matmul(rcplvs, np.asarray(split_line, dtype=float))
                frac = np.mod(frac, 1)
                self.trajectory.fractional_trajectory.append(frac.flatten())
            if split_line[0] in self.atom_list:
                self.trajectory.atom_name.append(split_line[0])
                name = True
                self.trajectory.total_atoms = self.trajectory.total_atoms + 1
                self.trajectory.atoms_in_history = self.trajectory.atoms_in_history + 1
        config.close()

    def _check_data(self):
        """
        Error handling function.
        """
        if self.trajectory.total_atoms == 0:
            raise ValueError("No Atoms of specified type exist within the CONFIG file")


class Archive():
    """
    The :py:class:`polypy.read.Trajectory` class evaluates the positions of all atoms in a ARCHIVE.

    Args:
        atom_list (:py:class:`list`): List of unique atom names in trajectory.
        datatype (:py:attr:`str`): Datatype of the original dataset e.g. DL_MONTE ARCHIVE.
    """
    def __init__(self, file, atom_list):
        self.file = file
        if os.path.isfile(self.file) is False:
            raise ValueError("File does not exist")
        self.atom_list = atom_list
        self.data_type = "DL_MONTE ARCHIVE"
        self.trajectory = Trajectory(self.atom_list,
                                     self.data_type)
        self.read_archive()
        self.trajectory._clean_data()

    def read_archive(self):
        """
        Read a DL_MONTE ARCHIVE file line by line and updates a :py:class:`polypy.read.Trajectory` object.
        """
        count = 0
        lv_count = 0
        timestep_atom_count = 0
        skipline = 1
        current_lv = []
        atom_name_encountered = False
        timestep = True
        archive = open(self.file, 'r')
        config_label = archive.readline().split()
        self.trajectory.timesteps = self.trajectory.timesteps + 1

        for line in archive:
            split_line = line.split()
            if lv_count == 3:
                current_lv = np.asarray(current_lv, dtype=float)
                rcplvs, length = ut.calculate_rcplvs(current_lv)
                self.trajectory.cell_lengths.append(length)
                self.trajectory.lv.append(current_lv)
                self.trajectory.reciprocal_lv.append(rcplvs)
                lv_count = 0
                skipline = 0
                timestep = False
                current_lv = []
            if lv_count < 3 and timestep is True and skipline == 2:
                current_lv.append(split_line)
                lv_count = lv_count + 1
            if atom_name_encountered:
                atom_name_encountered = False
                self.trajectory.cartesian_trajectory.append(split_line[:3])
                frac = np.matmul(rcplvs, np.asarray(split_line[:3], dtype=float))
                frac = np.mod(frac, 1)
                self.trajectory.fractional_trajectory.append(frac)
                timestep_atom_count = timestep_atom_count + 1
            if split_line[0] in self.atom_list:
                self.trajectory.atom_name.append(split_line[0])
                atom_name_encountered = True
                count = count + 1
            if split_line[0] == config_label[0]:
                self.trajectory.timesteps = self.trajectory.timesteps + 1
                timestep = True
                skipline = 1
                self.trajectory.atoms_at_timestep.append(timestep_atom_count)
                self.trajectory.record_number.append(self.trajectory.timesteps)
                self.trajectory.time.append(self.trajectory.timesteps)
                timestep_atom_count = 0
            elif timestep is True:
                skipline = 2

        archive.close()
