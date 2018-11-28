import numpy as np
from polypy import Write as wr
from polypy import Generic as ge
from polypy import Read as rd

class Density():
    '''
    The Density class will calculate the atomic density in one, two and both dimensions as well as within a specified plane. This
    class is ideal for observing how atomic density changes within a structure. 
    
    Parameters
    ----------
    data : dictionary containing the atomic trajectories, lattice vectors, timesteps and number of atoms. 
    '''

    def __init__(self, data, atom_type=None):
        self.data = data
        self.atom_type = atom_type
        if len(np.unique(self.data['label'])) > 1 and self.atom_type is None:
            print("WARNING: Multiple atom types detected - Splitting Coordinates")

        elif len(np.unique(self.data['label'])) > 1:
            self.data = rd.get_atom(self.data, self.atom_type)

    def one_dimensional_density(self, Bin=None, direction=None, output=None):
        '''
        one_dimensional_density - Calculate the atomic number density within one dimensional slices of a structure.
        Parameters
        ----------
        Bin             : Bin Value                       : Float                : Default - 0.1
        direction       : Direction normal to the bin     : String               : Default - x
        output          : Output file name                : String               : Default - 1D-Density.png
        
        Returns
        -------
        text file
        matplolib plot
        "Heatmap.png" - Contour plot - xy grid of total of atoms within each box
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
       
        c = self.data['trajectories'][:,val]
        b = (np.average(self.data['lv'][:,val]) / 2 ) 
        c = c + b
        x = ge.get_integer((np.amax(c)), Bin)
        bin_array = np.zeros((x))
        c.tolist()
    
        for j in range(0, self.data['trajectories'][:,val].size):
            
            plane = 0
            plane = ge.bin_choose(c[j], Bin)
            bin_array[plane] = bin_array[plane] + 1       
    
        x = np.arange( 0, ( bin_array.size ) )
        x = (x * Bin)  - b
        y = ( bin_array / self.data['timesteps'])
        wr.line_plot(x, y, "XCoordinate (" r'$\AA$' ")", "Number Density", output)
        return x, y
    
    def two_dimensional_density(self, box=None, direction=None, output=None):
        '''
        two_dimensional_density - 2D Atomic Density Analysis
        Parameters
        ----------
        box              : Box Value                             : Float                : Default : 0.1    
        direction        : Direction normal to the box           : String               : Default : x
        output           : output file name                      : String               : Default : 2D-Density.png
        log              : True for log plot, False for no log   : Boolean              : False
                    
        Returns
        -------
        matplolib plot
        "Heatmap.png" - Contour plot - xy grid of total of atoms within each box
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
        '''
        one_dimensional_density_sb - Calculate the total number of a given species within a 1D slice of a structure
        Parameters
        ----------
        ul           : Upper bin limit                         :    Float
        ll           : Lower bin limit                         :    Float
        direction    : Direction normal to slice               :    String          : Default - x
        
        Returns
        -------
        plane        :   Number of atomic species within bin   :    Integer
        '''
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
    
    def one_and_two_dimension_overlay(self, box=None, direction=None, output=None):
        '''
        one_and_two_dimension_overlay - Returns a plot of 1 dimensional density overlayed on two dimensional density
        Parameters
        ----------
        box        : Box Value                             : Float        : Default : 0.1    
        direction  : Direction normal to the box           : String       : Default : x
        output     : output file name                      : String       : Default : 2D-Density.png
        log        : True for log plot, False for no log   : Boolean      : False
                    
        Returns
        -------
        matplolib plot
        "Heatmap.png" - Contour plot - xy grid of total of atoms within each box
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