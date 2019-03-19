import numpy as np
from polypy import write as wr
from polypy import read as rd
from polypy import utils as ut


class Density():
    '''The Density class calculates the total number of atoms in one
    dimensional planes, two dimensional boxes and a combination of
    both. This class is ideal for observing how atomic density
    changes within a structure.

    Parameters
    ----------
    data : dict
        dictionary containing the atomic atom labels, trajectories,
        lattice vectors, number of timesteps and atoms.
    atom_type : list (optional)
        atoms to be analysed
    '''
    def __init__(self, data, atom_type=None):
        self.data = data
        self.atom_type = atom_type
        if len(np.unique(self.data['label'])) > 1 and self.atom_type is None:
            print("Multiple atom types detected - Splitting Coordinates")
        elif len(np.unique(self.data['label'])) > 1:
            self.data = rd.get_atom(self.data, self.atom_type)

    def one_dimensional_density(self, Bin=None, direction=None, output=None):
        '''Calculate the atomic number density within one dimensional
        slices of a structure.

        Parameters
        ----------
        Bin : float (optional)
            Size of one dimensional planes, perpendicular to a given direction.
        direction : string (optional)
            Direction perpendicular to the planes
        output : string  (optional)
            Name of the output graph

        Returns
        -------
        x : array like
            values corresponding to the bin coordinates.
        y : array like
            total number of species in each bin.
        '''
        if Bin is None:
            Bin = 0.1
        if direction is None:
            direction = "x"
        if output is None:
            output = "1D-Density.png"
        if direction == "x":
            val = 0
        elif direction == "y":
            val = 1
        elif direction == "z":
            val = 2
        c = self.data['trajectories'][:, val]
        b = (np.average(self.data['lv'][:, val]) / 2)
        c = c + b
        x = ut.get_integer((np.amax(c)), Bin)
        bin_array = np.zeros((x))
        c.tolist()

        for j in range(0, self.data['trajectories'][:, val].size):

            plane = 0
            plane = ut.bin_choose(c[j], Bin)
            bin_array[plane] = bin_array[plane] + 1

        x = np.arange(0, (bin_array.size))
        x = (x * Bin) - b

        return x, bin_array

    def two_dimensional_density(self, box=None, direction=None, output=None):
        '''Calculate the atomic number density within two dimensional
        boxes in a structure.

        Parameters
        ----------
        box : float (optional)
            Edge length of the boxes
        direction : string (optional)
            Direction pointing into the boxes
        output : string (optional)
            Name of the output plot

        Returns
        -------
        x : array like
            x axis coordinates
        y : array like
            y axis coordinates
        z : array like
            Total number of species in each box as a function of x and y.
        '''
        if box is None:
            box = 0.1
        if direction is None:
            direction = "x"
        if output is None:
            output = "2D-Density.png"
        if direction == "x":
            val = [1, 2]
        elif direction == "y":
            val = [0, 2]
        elif direction == "z":
            val = [0, 1]

        xc = self.data['trajectories'][:, val[0]]
        xc = xc + (np.average(self.data['lv'][:, [val[0]]]) / 2)
        yc = self.data['trajectories'][:, val[1]]
        yc = yc + (np.average(self.data['lv'][:, [val[1]]]) / 2)
        x = ut.get_integer(np.amax(xc), box)
        y = ut.get_integer(np.amax(yc), box)
        if x < y:
            x, y = y, x
            xc, yc = yc, xc
        bin_array = np.zeros(((y), (x)))
        xc = xc.tolist()
        yc = yc.tolist()

        for j in range(0, self.data['trajectories'][:, val[0]].size):

            xbox = 0
            ybox = 0
            xbox = ut.bin_choose(xc[j], box)
            ybox = ut.bin_choose(yc[j], box)
            bin_array[ybox, xbox] = bin_array[ybox, xbox] + 1

        bin_array = bin_array / self.data['timesteps']
        x = np.arange((x))
        y = np.arange((y))
        x = ((x * box)) - (np.average(self.data['lv'][:, [val[0]]]) / 2)
        y = ((y * box))
        z = bin_array + 0.001

        return x, y, z

    def one_dimensional_density_sb(self, ul=None, ll=None, direction=None):
        '''Calculate the total number of a given species within a
        1D slice of a structure.

        Parameters
        ----------
        ul : float
            Upper bin limit
        ll : float
            Lower bin limit
        direction : string
            direction perpendicular to the slice

        Returns
        -------
        plane : int
            Total number of species within specified bin.
        '''
        if direction is None:
            direction = "x"
        if direction == "x":
            val = 0
        elif direction == "y":
            val = 1
        elif direction == "z":
            val = 2

        c = self.data['trajectories'][:, val]
        c.tolist()
        plane = 0

        if ul is None:
            ul = np.amax(c)
        if ll is None:
            ll = np.amin(c)

        for j in range(0, self.data['trajectories'][:, val].size):

            if c[j] > ll and c[j] < ul:
                plane = plane + 1

        return plane

    def one_and_two_dimension_overlay(self, box=None, direction=None,
                                      output=None):
        '''Combination of the one and two dimensional density functions.
        Calculates both total number of atoms in one and two dimensions
        and overlays them in one plot.

        Parameters
        ----------
        box : float (optional)
            size of the boxes and bins
        direction : string (optional)
            direction perpendicular to the planes
        output : string (optional)
            Output file name.

        Returns
        -------
        x : array like
            X axis coordinates.
        y : array like
            Y axis coordinates.
        od_array : array like
            Total number of species within each plane.
        td_array : array like
            Total number of species within each box.
        '''
        if box is None:
            box = 0.1
        if direction is None:
            direction = "x"
        if output is None:
            output = "Combined-Density.png"
        if direction == "x":
            val = [1, 2]
        elif direction == "y":
            val = [0, 2]
        elif direction == "z":
            val = [0, 1]

        xc = self.data['trajectories'][:, val[0]]
        xc = xc + (np.average(self.data['lv'][:, [val[0]]]) / 2)
        yc = self.data['trajectories'][:, val[1]]
        yc = yc + (np.average(self.data['lv'][:, [val[1]]]) / 2)
        x = ut.get_integer(np.amax(xc), box)
        y = ut.get_integer(np.amax(yc), box)
        od_array = np.zeros((x))
        td_array = np.zeros(((y), (x)))
        xc = xc.tolist()
        yc = yc.tolist()

        for j in range(0, self.data['trajectories'][:, val[0]].size):

            xbox = 0
            ybox = 0
            xbox = ut.bin_choose(xc[j], box)
            ybox = ut.bin_choose(yc[j], box)
            od_array[xbox] = od_array[xbox] + 1
            td_array[ybox, xbox] = td_array[ybox, xbox] + 1

        td_array = td_array / self.data['timesteps']
        od_array = od_array / self.data['timesteps']
        x = np.arange((x))
        y = np.arange((y))
        x = ((x * box)) - (np.average(self.data['lv'][:, [val[0]]]) / 2)
        y = ((y * box))
        td_array = td_array + 0.001
        od_array = od_array + 0.001

        return x, y, td_array, od_array
