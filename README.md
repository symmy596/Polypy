# PolyPy 

A short program to analyse the outputs from DL_POLY MD simultations.

Author : Adam Symington   
Date   : 21/06/2018  


When this was committed only I and God himself understood how it worked. 


### Basic usage :

As a first step you need to read in the trajectory that you want for a given atom. This will 
return the number of atoms, number of timesteps and the trajectories. 

Natoms, NConfigs, Coords = rd.ReadHistory("Filename", "Atom")  
  
Filename - Name of History file  
Atom - Atom of interest  

Natoms - Number of Atoms  
NConfigs - Number of Timesteps
Coords - Trajectories  

Once the coordinates are read in then they can be fed into the various other bits of functionality.  

### Functionality

####1) 1 dimensional density plot  

ta.one_dimensional_density(Coords, NAtoms, NConfigs, Vec, Bin, "x")

<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/1D-Density.png">
</p>

####2) 2 dimensional density plot

ta.two_dimensional_density(Coords, NAtoms, NConfigs, Vec, Box, 'z')

<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/2D-Density.png">
</p>


####3) Single MSD run

ta.msd(Coords, Vec, NConfigs, NAtoms)


####4) MSD within a specific region of the system
- This needs work

ta.plane_msd(Coords, NConfigs, NAtoms, UL, LL, Vec)


####5) Smoothed MSD

ta.smooth_msd(Coords, Vec, Runs, NConfigs, NAtoms)


####6) MSD that plots diffusion coefficient of each atom against its average position 

ta.pmsd(Coords, Vec, NConfigs, NAtoms, Bin)

####7) System Volume

volume, time = ta.system_volume(lv, NConfigs, timestep)


