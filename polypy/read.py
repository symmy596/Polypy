import os as os
import sys as sys
import numpy as np
from polypy import utils as ut

class Trajectory():

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

    def _list_to_arrays(self):
        self.cartesian_trajectory = np.asarray(self.cartesian_trajectory, dtype=float)
        self.fractional_trajectory = np.asarray(self.fractional_trajectory, dtype=float)
        self.atom_name = np.asarray(self.atom_name, dtype=str)
        self.lv = np.asarray(self.lv, dtype=float)
        self.reciprocal_lv = np.asarray(self.reciprocal_lv, dtype=float)
        self.cell_lengths = np.asarray(self.cell_lengths, dtype=float)

    def get_config(self, timestep):
        if self.data_type == "DL_POLY CONFIG":
            raise ValueError("Only one timestep was found")
        
        config_trajectory = Trajectory(self.atom_list, "DL_POLY CONFIG")
        config_trajectory.cartesian_trajectory = np.split(self.cartesian_trajectory, self.timesteps)[timestep]
        config_trajectory.fractional_trajectory = np.split(self.fractional_trajectory, self.timesteps)[timestep]
        config_trajectory.atom_name = np.split(self.atom_name, self.timesteps)[timestep]
        config_trajectory.reciprocal_lv = self.reciprocal_lv[timestep]
        config_trajectory.lv = self.lv[timestep]
        config_trajectory.cell_lengths = np.split(self.cell_lengths, self.timesteps)[timestep]
        config_trajectory.atoms_in_history = self.total_atoms
        config_trajectory.timesteps = 1
        config_trajectory.total_atoms = self.total_atoms
        return config_trajectory

    def get_atom(self, atom):
        if np.unique(self.atom_list).size == 1:
            raise ValueError("Only one atom type was found")

        atom_trajectory = Trajectory([atom], self.data_type)

        for i in range(0, self.atom_name.size):
            if self.atom_name[i] == atom:
                atom_trajectory.cartesian_trajectory.append(self.cartesian_trajectory[i])
                atom_trajectory.fractional_trajectory.append(self.fractional_trajectory[i])
                atom_trajectory.atom_name.append(self.atom_name[i])
                atom_trajectory.atoms_in_history = atom_trajectory.atoms_in_history + 1

        atom_trajectory.reciprocal_lv = self.reciprocal_lv
        atom_trajectory.lv = self.lv
        atom_trajectory.cell_lengths = self.cell_lengths
        atom_trajectory.timesteps = self.timesteps
        atom_trajectory.total_atoms = atom_trajectory.atoms_in_history / self.timesteps
        atom_trajectory._list_to_arrays()

        return atom_trajectory

class History():

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
        self.trajectory._list_to_arrays()

    def read_history(self):
        '''Read a DL_POLY HISTORY file.
        '''
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
        history.close()

    def _check_data(self):
        if self.trajectory.total_atoms == int(self.trajectory.total_atoms) is False:
            raise ValueError("The total number of atoms is not constant across each timestep")
        if self.trajectory.total_atoms == 0:
            raise ValueError("No Atoms of specified type exist within the HISTORY file")

class Config():

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
        self.trajectory._list_to_arrays()


    def read_config(self):

        config = open(self.file, 'r')
        name = False

        title = config.readline()
        stuff = config.readline()
        lv = []
        for i in range(0, 3):
            l = config.readline()
            lv.append(l.split())
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
                frac = np.mod(frac, 1 )
                self.trajectory.fractional_trajectory.append(frac.flatten())
            if split_line[0] in self.atom_list:
                self.trajectory.atom_name.append(split_line[0])
                name = True
                self.trajectory.total_atoms = self.trajectory.total_atoms + 1
                self.trajectory.atoms_in_history = self.trajectory.atoms_in_history + 1

        config.close()

    def _check_data(self):
        if self.trajectory.total_atoms == 0:
            raise ValueError("No Atoms of specified type exist within the CONFIG file")
    
def read_archive(file, atom_list):
    '''Read a DL_MONTE ARCHIVE file
    Parameters
    ----------
    file : str
        Name of the dlmonte ARCHIVE file
    atom_list : list
        list of atoms types to be read

    Returns
    -------
    data : dict
        Dictionary containing the atom labels, trajectories,
        lattice vectors, timesteps and number of atoms
    '''
    if os.path.isfile(file):
        trajectories = []
        atname = []
        lv = []
        natoms_at_timestep = []
        count = 0
        c = 0
        atom_count = 0
        timesteps = 1
        skipline = 1
        name = False
        tstep = True
        archive = open(file, 'r')
        config_label = archive.readline()
        config_label = config_label.split()
        for line in archive:
            x = line.split()
            if c == 3:
                c = 0
                skipline = 0
                tstep = False
            if c < 3 and tstep is True and skipline == 2:
                lv.append(line.split())
                c = c + 1
            if name:
                name = False
                trajectories.append(line.split())
                atom_count = atom_count + 1
            if x[0] in atom_list:
                atname.append(x[0])
                name = True
                count = count + 1
            if x[0] == config_label[0]:
                timesteps = timesteps + 1
                tstep = True
                skipline = 1
                natoms_at_timestep.append(atom_count)
                atom_count = 0
            elif tstep is True:
                skipline = 2

        trajectories = np.asarray(trajectories, dtype=float)
        atname = np.asarray(atname, dtype=str)
        lv = np.asarray(lv, dtype=float)
        natoms = count / timesteps
        natoms = int(natoms)
        vec = np.array([])
        lv = np.split(lv, timesteps)

        for i in range(0, timesteps):

            vec = np.append(vec, (lv[i].sum(axis=0)))

        lv = np.reshape(vec, (timesteps, 3))
        data = {'label': atname,
                'trajectories': trajectories,
                'lv': lv,
                'timesteps': timesteps,
                'natoms': natoms,
                "atoms_per_timestep": natoms_at_timestep}

    else:
        print("File cannot be found")
        sys.exit(0)

    if natoms == 0:
        print("No Atoms of specified type exist within the ARCHIVE file")
        sys.exit(0)

    archive.close()
    return data




def get_atom(data, atom):
    """Reads through the dictionary returned from the reading functions
    and returns a given atom.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    atom : str
        atom type to be found

    Returns
    -------
    atom_data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    """
    if len(np.unique(data['label'])) == 1:
        return data
    frac_coords = []
    coords = []
    count = 0
    for i in range(0, data['label'].size):
        if data['label'][i] == atom:
            coords.append(data['trajectories'][i])
            frac_coords.append(data['frac_trajectories'][i])
            count = count + 1

    coords = np.asarray(coords)
    frac_coords = np.asarray(frac_coords)
    natoms = int(count / data['timesteps'])
    atom_data = {'label': atom,
                 'trajectories': coords,
                 'frac_trajectories': frac_coords,
                 'lv': data['lv'],
                 'lengths': data['lengths'],
                 'timesteps': data['timesteps'],
                 'natoms': natoms}

    return atom_data


def get_config(data, timestep):
    """Returns single config from entire trajectory

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : int
        Timestep of desired config

    Returns
    -------
    array like
        Coordinates of desired config
    """
    configs = np.split(data['trajectories'], data['timesteps'])
    return configs[timestep]

def get_trajectory(data, atom):
    configs = np.split(data['trajectories'], data['timesteps'])
    configs = np.asarray(configs)
    return configs[:,atom]
