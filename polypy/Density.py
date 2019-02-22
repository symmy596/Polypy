import numpy as np
from polypy import Write as wr
from polypy import Generic as ge
from polypy import Read as rd

class Density():
    """Density system.
    This class is designed to take atomic trajectories and calculate the
    average number density of that species within bins in a given
    direction (One dimension) or within boxes in two given directions
    (Two dimensions). Both the one dimensional and two dimensional
    can be calculated and displayed as one.
    
    Parameters
    ----------
    data : dictionary 
        Dictionary containing the atomic trajectories, lattice vectors,
        timesteps and number of atoms. 
    atom_type : optional
        Specify the atom to retrieve the density for.
    """

    def __init__(self, data, atom_type=None):
        self.data = data
        self.atom_type = atom_type
        if len(np.unique(self.data['label'])) > 1 and self.atom_type is None:
            print("Multiple atom types detected - Splitting Coordinates")
        elif len(np.unique(self.data['label'])) > 1:
            self.data = rd.get_atom(self.data, self.atom_type)

    def one_dimensional_density(self, bin=None, direction=None, output=None):
        """Calculate the number density within bins perpendicular 
        to a given direction within a structure.

        Parameters
        ----------
        bin : float (optional)
            bin size in a given direction.
        direction : string (optional)
            direction perpendicular to the bins.
        output : string (optional)
            name of output plot.
        
        Returns
        -------
        x : float
            bin coordinates
        y : float
            number denisty in each bin
        """
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
        c = self.data['trajectories'][:,val]
        b = (np.average(self.data['lv'][:,val]) / 2 ) 
        c = c + b
        x = ge.get_integer((np.amax(self.data['lv'][:,val])), Bin)
        bin_array = np.zeros((x))
        c.tolist()
        for j in range(0, self.data['trajectories'][:,val].size):  
            plane = 0
            plane = ge.bin_choose(c[j], Bin)
            bin_array[plane] = bin_array[plane] + 1       
        x = np.arange( 0, ( bin_array.size ) )
        x = (x * Bin)  - b
        y = ( bin_array / self.data['timesteps'])
        wr.line_plot(x, y, "XCoordinate (" r'$\AA$' ")", "Number Density", 
                     output)
        return x, y
    
    def two_dimensional_density(self, box=None, direction=None, output=None):
        """Calculate the number density within boxes.

        Parameters
        ----------
        box : float (optional)
            box size.
        direction : string (optional)
            direction perpendcular to boxes.
        output : string (optional)
            name of output plot.
        
        Returns
        -------
        x : float
            box coordinates in x
        y : float
            box coordinates in y
        z : float
            number denisty in box
        """
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
        xc = self.data['trajectories'][:,val[0]]
        xc = xc + ( np.average(self.data['lv'][:,[val[0]]]) / 2 )             
        yc = self.data['trajectories'][:,val[1]]
        yc = yc + ( np.average(self.data['lv'][:,[val[1]]]) / 2 ) 
        x = ge.get_integer(np.amax(xc), box)
        y = ge.get_integer(np.amax(yc), box)
        bin_array = np.zeros(((y), (x)))
        xc = xc.tolist()
        yc = yc.tolist()
        for j in range(0, self.data['trajectories'][:,val[0]].size):        
            xbox = 0
            ybox = 0
            xbox = ge.bin_choose(xc[j], box)
            ybox = ge.bin_choose(yc[j], box)
            bin_array[ybox, xbox] = bin_array[ybox, xbox] + 1       
        bin_array = bin_array / self.data['timesteps']
        x = np.arange((x))
        y = np.arange((y))
        x = ((x * box)) - (np.average(self.data['lv'][:,[val[0]]]) / 2 )
        y = ((y * box))
        z = bin_array + 0.001
        wr.contour_plot(x, y, z, output)
        return x, y, z

    def one_dimensional_density_sb(self, ul=None, ll=None, direction=None):
        """Calculate the number density within a single bin perpendicular
        to a given direction within a structure.

        Parameters
        ----------
        ul : float (optional)
            upper edge of bin.
        ll : float (optional)
            lower edge of bin.
        direction : string (optional)
            direction perpendicular to the bin.            
        
        Returns
        -------
        plane : float
            Total number of species in bin
        """
        if direction is None:
            direction = "x" 
        if direction == "x":
            val = 0
        elif direction == "y":
            val = 1
        elif direction == "z":
            val = 2 
        c = self.data['trajectories'][:,val]
        c.tolist()
        plane = 0
        if ul is None:
            ul = np.amax(c)
        if ll is None:
            ll = np.amin(c)
        for j in range(0, self.data['trajectories'][:,val].size):      
            if c[j] > ll and c[j] < ul:
                plane = plane + 1       
        filename = "1D-Density-" + (str(ul)) + " - " + (str(ll))                
        wr.one_dimensional_density_sb_output(plane, ul, ll, filename)   
        return plane 
    
    def one_and_two_dimension_overlay(self, box=None, direction=None,
                                      output=None):
        """Calculate the number density within boxes within a structure. A
        one dimensional version is displayed over the 2D plot

        Parameters
        ----------
        bin : float (optional)
            bin size in a given direction.
        direction : string (optional)
            direction perpendicular to the bins.
        output : string (optional)
            name of output plot.
        
        Returns
        -------
        x : float
            bin coordinates
        y : float
            number denisty in each bin
        """
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
    
        xc = self.data['trajectories'][:,val[0]]
        xc = xc + ( np.average(self.data['lv'][:,[val[0]]]) / 2 )             
        yc = self.data['trajectories'][:,val[1]]
        yc = yc + ( np.average(self.data['lv'][:,[val[1]]]) / 2 ) 
        x = ge.get_integer(np.amax(xc), box)
        y = ge.get_integer(np.amax(yc), box)
        od_array = np.zeros((x))
        td_array = np.zeros(((y), (x)))
        xc = xc.tolist()
        yc = yc.tolist()
        for j in range(0, self.data['trajectories'][:,val[0]].size):        
            xbox = 0
            ybox = 0
            xbox = ge.bin_choose(xc[j], box)
            ybox = ge.bin_choose(yc[j], box)
            od_array[xbox] = od_array[xbox] + 1
            td_array[ybox, xbox] = td_array[ybox, xbox] + 1       
        td_array = td_array / self.data['timesteps']
        od_array = od_array / self.data['timesteps']
        x = np.arange((x))
        y = np.arange((y))
        x = ((x * box)) - (np.average(self.data['lv'][:,[val[0]]]) / 2 )
        y = ((y * box))
        td_array = td_array + 0.001
        od_array = od_array + 0.001
        wr.combined_density_plot(x, y, od_array, td_array, output)

def charge_density(densities, charges, bin_volume):
    """Calculates the charge density from given number densities
    and atomic charges.

    Parameters
    ----------
    densities : list
        list of numpy arrays containing atomic densities.
    charges : list
        list of charges corresponding to species in densities.
    bin_volume : float
        Volume of bins

    Returns
    -------
    charge_density : array like
        charge density
    """
    pass
#    charge_density = densities * charges
    

def electric_field(charge_density, bin_positions):
    """Calculates the electric field from the charge density.

    Parameters
    ----------
    charge_density : array like
        array of charge densities
    bin_positions : array like
        array of bins

    Returns
    -------
    electric_field : array like
        array of electric field
    """
    pass

def electrostatic_potential(electric_field, bin_positions):
    """Calculates the electrostatic potential from the charge density.

    Parameters
    ----------
    electric_field: array like
        array of electric field
    bin_positions : array like
        array of bins

    Returns
    -------
    electrostatic_potential : array like
        array of electrostatic potentials
    """
    pass

