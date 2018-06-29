# PolyPy 

A short program to analyse the outputs from DL_POLY MD simultations.

Author : Adam Symington   
Date   : 21/06/2018  
  
This program should in theory be able to calculate atomic density in 1 and 2D, do various mean squared displacement calculations and do a volume analysis.  
  

### Basic usage :

As a first step you need to read in the trajectory that you want for a given atom. This will 
return the number of atoms, number of timesteps and the trajectories. 

natoms, timesteps, trajectories, lv = rd.ReadHistory("Filename", "Atom")  
  
Filename - Name of History file  
Atom - Atom of interest  

natoms - Number of Atoms  
timesteps - Number of Timesteps
trajectories - Trajectories  
lv = Lattice vectors at each timestep

Once the coordinates are read in then they can be fed into the various other bits of functionality.  

### Functionality

#### 1) One Dimensional Density Plot  

This function does a simple one dimensional density calculation for a given atom. 

ta.one_dimensional_density(trajectories, natoms, timesteps, lv, bin, "x", output)

trajectories - Trajectories   
natoms - Number of Atoms  
timesteps - Number of timesteps  
lv - Lattice Vectors  
bin - Bin size - default = 0.1  
"x" - Direction - default = "x"  
output - output file name - default = 1D-Density.png  

  
Returns a plot showing atomic number density relative to the position within the cell.   
<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/1D-Density.png">
</p>

#### 2) Two dimensional density plot

This function does a simple two dimensional density calculation for a given atom.

ta.two_dimensional_density(trajectories, natoms, timesteps, lv, box, 'z', output)

trajectories - Trajectories
natoms - Number of Atoms
timesteps - Number of timesteps
lv - Lattice Vectors
box - Box size - default = 0.1
"x" - Direction normal to the boxes - default = "x"  
output - output file name - default = 2D-Density.png   
  
Returns a heatmap of atomic positions in 2D. 
  
<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/2D-Density.png">
</p>


#### 3) Single MSD run

ta.msd(trajectories, timesteps, natoms, timestep, lv)

trajectories - Trajectories  
natoms - Number of Atoms  
timesteps - Number of timesteps  
lv - Lattice Vectors  

[Link to my blog post on MSD](http://people.bath.ac.uk/ars44/blog_posts/post_1.html)

<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/MSD.png">
</p>

#### 4) System Volume

Simple function to read the lattice vectors at each timestep and calculate the cell volume.  

volume, time = ta.system_volume(lv, timesteps, timestep)
  

  
<p align="center">
  <img width="460" height="300" src="https://github.com/symmy596/PolyPy/blob/master/Plots/Volume.png">
</p>


#### 5) MSD within a specific region of the system


ta.plane_msd(trajectories, timesteps, natoms, UL, LL, direction, lv, timestep)


#### 6) Smoothed MSD

ta.smooth_msd(trajectories, runs, timestep, natoms, lv, timestep)


#### 7) MSD that plots diffusion coefficient of each atom against its average position 

ta.pmsd(trajectories, lv, timesteps, natoms, bin, direction)



