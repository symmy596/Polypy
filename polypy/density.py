import numpy as np
from polypy import read as rd
from polypy import utils as ut


class Density():
    '''The Density class calculates the total number of atoms in one
    dimensional planes, two dimensional boxes and a combination of
    both. This class is ideal for observing how particle density
    changes within a structure.

    Parameters
    ----------
    data : obj
        Class object containing the trajectory
    atom_type : list (optional)
        atoms to be analysed
    '''
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

        self.x_lim = np.ceil(np.amax(self.lengths[:, 0]) / self.histogram_size).astype(int)
        self.y_lim = np.ceil(np.amax(self.lengths[:, 1]) / self.histogram_size).astype(int)
        self.z_lim = np.ceil(np.amax(self.lengths[:, 2]) / self.histogram_size).astype(int)
        self.x = (np.arange(0, self.x_lim)) * self.histogram_size
        self.y = (np.arange(0, self.y_lim)) * self.histogram_size
        self.z = (np.arange(0, self.z_lim)) * self.histogram_size

    def update_map(self, position):

            xbox = (position[0] * self.x_lim).astype(int)
            ybox = (position[1] * self.y_lim).astype(int)
            zbox = (position[2] * self.z_lim).astype(int)
            self.coords_map[ybox, xbox, zbox] = self.coords_map[ybox, xbox, zbox] + 1

    def build_map(self):

        positions = self.data.fractional_trajectory.tolist()

        for position in positions:
            
            self.update_map(position)

    def one_dimensional_density(self, direction="x"):
        '''Calculate the particle density within one dimensional
        slices of a structure.

        Parameters
        ----------
        histogram_width : float (optional)
            Width of histograms, perpendicular to a given direction.
        direction : string (optional)
            Direction perpendicular to the histograms

        Returns
        -------
        x : array like
            values corresponding to the histograms coordinates.
        histograms : array like
            total number of species in each histograms.
        '''
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

    def two_dimensional_density(self, box=0.1, direction="x"):
        '''Calculate the atomic number density within two dimensional
        boxes in a structure.

        Parameters
        ----------
        box : float (optional)
            Edge length of the boxes
        direction : string (optional)
            Direction pointing into the boxes

        Returns
        -------
        x : array like
            x axis coordinates
        y : array like
            y axis coordinates
        z : array like
            Total number of species in each box as a function of x and y.
        '''
        if direction == "x":
            val = 0
            x = self.y
            y = self.z
        elif direction == "y":
            val = 1
            x = self.x
            y = self.z
        elif direction == "z":
            val = 2
            x = self.x
            y = self.y
        z = np.sum(self.coords_map, axis=val)
        box_volume = 0.1 * 0.1 * np.mean(self.data.cell_lengths[:,val])
        return x, y, z, box_volume


    def one_dimensional_density_sb(self, ul, ll, direction="x"):
        '''Calculate the total number of a given species within a
        1D slice of a structure.

        Parameters
        ----------
        ul : float
            Upper bin limit
        ll : float
            Lower bin limit
        direction : string (optional)
            direction perpendicular to the slice

        Returns
        -------
        plane : int
            Total number of species within specified bin.
        '''
        if direction == "x":
            val = 0
        elif direction == "y":
            val = 1
        elif direction == "z":
            val = 2

        c = self.coords[:, val]
        c.tolist()
        plane = 0

        for j in range(0, self.coords[:, val].size):

            if c[j] > ll and c[j] < ul:
                plane = plane + 1

        return plane

def regional_residence_time(data, ul, ll, direction, timestep):
    if direction == "x":
        val = 0
    elif direction == "y":
        val = 1
    elif direction == "z":
        val = 2

    times = np.array([])
    for i in range(data['natoms']):
        trajectory = rd.get_trajectory(data, i)
        in_region = False
        time_in_box = 0
        previously_in = False
        for j in range(trajectory[:,0].size):
           # print(trajectory[i,val])

            if trajectory[j,val] > ll and trajectory[j,val] < ul:
                if previously_in == True:  
                    time_in_box = time_in_box + 1

                if in_region == False:
                    previously_in = True

                in_region = True            
            elif trajectory[j,val] < ll or trajectory[j,val] > ul:
                if in_region == True:
                    in_region = False
                    previously_in = False
                    if time_in_box > 10:
                        times = np.append(times, time_in_box)
                    time_in_box = 0
            
           # print(trajectory[j,val], in_region, previously_in, j, i, time_in_box)
    print(times)
    steps = np.average(times)
    times = steps * timestep
    return steps, times