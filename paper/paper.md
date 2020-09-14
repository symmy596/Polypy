---
title: 'polypy - '
tags:
- Chemistry
- Physics
- Materials Science
- Solid State Chemistry
- Simulation
- Molecular Dynamics
- Monte Carlo
authors:
- name: Adam R. Symington
  orcid: 0000-0001-6059-497X
  affiliation: "1"
affiliations:
- name: Department of Chemistry, University of Bath
  index: 1
date: 11 September 2020
bibliography: paper.bib
---

# Summary

A large number of research questions in solid state chemistry can be addressed using molecular dynamics and Monte Carlo simulations. These simulations allow many material properties to be calculated for direct comparison with experiment. These include, the diffusion coefficients and ionic conductivities, charge density, electric field and electrostatic potential. For example, the diffusion coefficient and ionic conductivity are of particular importance for the study of battery materials (e.g. Li-ion / Na-ion diffusion), solid oxide fuel cell materials (e.g. O-ion diffusion) and many other applications in solid state chemistry. The charge density, electric field and electrostatic potential are of interest to problems relating to interfaces in solid state chemistry, e.g. grain boundaries and surfaces. Finally, calculating the distribution of defects in a material is useful for the study of segregation behaviour.

A molecular dynamics trajectory is a snapshot of the positions of each atom as a function time e.g. the trajectory of a single atom would show, sequentially, all of the postiions occupied by that atom throughout the simulaton. A Monte Carlo trajectory is similar although the simulation is not time resolved and the atom positions are simply a function of simulation step, not simulation timestep. The positions of the atoms allows the number density of each atom to be calculated and from these, the electrostatic potential ($\Delta_{\psi}(z)$), electric field ($E(z)$) and charge density ($\rho_q(z)$) can be calculated according to

\begin{align}
\Delta_{\psi}(z) = \int_{z_{0}}^{z} E(z')dz',
\end{align}

where, 

\begin{align}
E(z) = \frac{1}{- \epsilon_{0}} \int_{z_{0}}^{z} \rho_{q}(z')dz',
\end{align}

where, 

\begin{align}
\rho_q(z) = \sum_{i} q_i \rho_i(z),
\end{align}

$\rho_{i}$ is the density of atom i and $q_{i}$ is its charge.  

The atomic trajectories as a function of time allow diffusion coefficients to be calculated using a mean squared displacement. Using a dye molecule in water as an example, the motion of a dye molecule is not simple. As it moves it is jostled by collisions with other molecules, preventing it from moving in a straight path. If the path is examined in close detail, it will be seen to be a good approximation to a random walk. In mathematics a random walk is a series of steps, each taken in a random direction. According to Einstein the mean square of the distance travelled by a particle following a random walk is proportional to the time elapsed.

\begin{align}
\Big \langle r_{i}^{2} \big \rangle & = 6 D_t + C 
\end{align}

where

\begin{align}
\Big \langle r_{i}^{2} \big \rangle = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big>.
\end{align}

where $\Big \langle r^2 \big \rangle$ is the mean squared distance, t is time, $D_t$ is the diffusion rate and C is a constant. If $\Big \langle r_{i}^{2} \big \rangle$ is plotted as a function of time, the gradient of the curve obtained is equal to 6 times the self-diffusion coefficient of particle i.

The `polypy` code is designed to solve the following problems.

- Read DL_POLY and DL_MONTE trajectories.
- Calculate the number density of all species in a trajectory in one and two dimensions.
- Calculate the charge density in one and two dimensions.
- Calculate the electric field and electrostatic potential in one dimension.
- Calculate the mean squared displacement for a given atom and use this to calculate the diffusion coefficient and ionic conductivity.
- Calculate the volume as a function of simulation timestep.
- Generate publication ready plots.

`polypy` has been used to study Li-ion diffusion in lithium lanthanum titanate [@LLTO], oxygen diffusion in cerium oxide [@CeO2] and oxygen diffusion / cation migration in uranium oxide [@UO2]. 

<p align="center"> 
<img src="fig_1.png" width="100%"/>
</p>

# `polypy`

`polypy` is a Python module for analysing molecular dynamics and Monte Carlo trajectories generated from the DL_POLY and DL_MONTE codes.
The code reads DL_POLY HISTORY and CONFIG files, and DL_MONTE ARCHIVE files and stores the data in a `polypy.read.Trajectory` object which is then used by the various data analysis modules.
The `polypy.density.Density` class generates a three dimensional grid and counts the number of times a given atom spends in each grid point during the simulation. This is then used to generate the number density of a given atom in one and two dimension. From here, the charge density in one and two dimensions, the electric field in one dimension and electrostatic potential in one dimension can be calculated using the `polypy.analysis` module. 
The `polypy.msd` module performs an mean squared displacement calculation. From the mean squared displacement, the three, two and one dimensional diffusion coefficient and ionic coefficient can be calculated. A module allowing easy generation of publication plots from the calculated data is available. The outputs are returned in a sensible form, allowing further manipulation and plotting.
`polypy` is aimed towards theoretical solid state physicist who have a basic familiarity with Python.
The repository contains examples of the core functionality as well as tutorials, implemented in Jupyter notebooks to explain the full theory. Furthermore, a detailed description of theory is also available within the documentation.

# Acknowledgements
  
ARS would like to thank Andrew R. McCluskey, Benjamin Morgan and Stephen C. Parker for their help and guidance. This package was written during a PhD funded by AWE and EPSRC (EP/R010366/1). 

# References