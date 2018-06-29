# PolyPy 

A short program to analyse the outputs from DL_POLY MD simultations.

Author : Adam Symington   
Date   : 21/06/2018  


When this was committed only I and God himself understood how it worked. 


### Basic usage :

As a first step you need to read in the trajectory that you want for a given atom. This will 
return the number of atoms, number of timesteps and the trajectories. 

Natoms, NConfigs, Coords, lv = rd.ReadHistory("Filename", "Atom")  
  
Filename - Name of History file  
Atom - Atom of interest  

Natoms - Number of Atoms  
NConfigs - Number of Timesteps
Coords - Trajectories  
lv = Lattice vectors at each timestep

Once the coordinates are read in then they can be fed into the various other bits of functionality.  

### Functionality

#### 1) One Dimensional Density Plot  

This function does a simple one dimensional density calculation for a given atom. 

ta.one_dimensional_density(Coords, NAtoms, NConfigs, lv, Bin, "x", output)

Coords - Trajectories   
NAtoms - Number of Atoms  
NConfigs - Number of timesteps  
lv - Lattice Vectors  
Bin - Bin size - default = 0.1  
"x" = Direction - default = "x"  
output - output file name - default = 1D-Density.png  

  
Returns a plot showing atomic number density relative to the position within the cell.   
<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/1D-Density.png">
</p>

#### 2) Two dimensional density plot

This function does a simple two dimensional density calculation for a given atom.

ta.two_dimensional_density(Coords, NAtoms, NConfigs, lv, Box, 'z', output)

Coords - Trajectories  
NAtoms - Number of Atoms  
NConfigs - Number of timesteps  
lv - Lattice Vectors  
Box - Box size - default = 0.1  
"x" - Direction normal to the boxes - default = "x"  
output - output file name - default = 2D-Density.png   
  
Returns a heatmap of atomic positions in 2D. 
  
<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/2D-Density.png">
</p>


#### 3) Single MSD run

ta.msd(Coords, lv, NConfigs, NAtoms)

Coords - Trajectories  
NAtoms - Number of Atoms  
NConfigs - Number of timesteps  
lv - Lattice Vectors  

[Link to my blog post on MSD](http://people.bath.ac.uk/ars44/blog_posts/post_1.html)

<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/MSD.png">
</p>

#### 7) System Volume

Simple function to read the lattice vectors at each timestep and calculate the cell volume.  

volume, time = ta.system_volume(lv, NConfigs, timestep)
  

  
<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/Volume.png">
</p>


#### 5) MSD within a specific region of the system


ta.plane_msd(Coords, NConfigs, NAtoms, UL, LL, Vec)


#### 6) Smoothed MSD

ta.smooth_msd(Coords, Vec, Runs, NConfigs, NAtoms)


#### 7) MSD that plots diffusion coefficient of each atom against its average position 

ta.pmsd(Coords, Vec, NConfigs, NAtoms, Bin)



