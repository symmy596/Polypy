---
title: 'polypy - Analysis Tools for Solid State Molecular Dynamics and Monte Carlo Trajectories'
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
date: 29 September 2020
bibliography: paper.bib
---

# Summary

A large number of research questions in solid state chemistry can be addressed using molecular dynamics and Monte Carlo simulations. These simulations allow many material properties to be calculated for direct comparison with experiment. These include, the diffusion coefficients, ionic conductivities, charge density, electric field and electrostatic potential. The diffusion coefficient and ionic conductivity are of particular importance for the study of battery materials (e.g. Li-ion / Na-ion diffusion [@LLTO]), solid oxide fuel cell materials (e.g. O-ion diffusion [@CeO2]) and many other applications in solid state chemistry. The charge density, electric field and electrostatic potential are of interest to problems relating to interfaces in solid state chemistry, e.g. space charge theory .[@Maier_BerBunsenges1984; @Maier_JPhysChemSol1985; @Maier_SolStatIonics2003; @ChiangEtAl_ApplPhysLett1996; @KimAndMaier_JElectrochemSoc2002] Finally, calculating the distribution of defects in a material is useful for the study of segregation behaviour [@UO2; @CeO2] or adsorption behaviour.[@Nora]

A molecular dynamics trajectory is a snapshot of the positions occupied by each atom in the simulation as a function time. For example, the trajectory of a single atom would show, sequentially, all of the postiions occupied by that atom throughout the simulaton. A Monte Carlo trajectory is similar although the simulation is not time resolved and the atom positions are simply a function of simulation step, not simulation timestep. The positions of the atoms allow the particle density of each atom to be determined and from these, the electrostatic potential, electric field and charge density can be calculated. A mean squared displacement can be performed on the molecular dynamics trajectories and from these, the diffusion coefficients and ionic conductivites can be calculated. Diffusion coefficients and ionic conductivites can then be used to estimate the activation energy for diffusion using the Arrhenius relationship. 

The `polypy` code is designed to solve the following problems.

- Read DL_POLY [@Smith] and DL_MONTE [@Purton] trajectories.
- Calculate the particle density of all species in a trajectory in one and two dimensions.
- Calculate the charge density in one and two dimensions.
- Calculate the electric field and electrostatic potential in one dimension.
- Calculate the mean squared displacement for a given atom and use this to calculate the diffusion coefficient and ionic conductivity.
- Calculate the volume as a function of simulation timestep.
- Generate publication ready figures.

`polypy` has been used to study Li-ion diffusion in lithium lanthanum titanate, [@LLTO] oxygen diffusion and cation migration in both uranium oxide and cerium oxide, [@UO2; @CeO2] thus there is a clear research application. 

<p align="center"> 
<img src="fig_1.png" width="100%"/>
<figcaption>Figure 1 - (a) Mean squared displacement for fluorine diffusion in calcium fluoride. The msd has been plotted in one, two and three dimensions. (b) The evolution of the system volume during a molecular dynamics simulation. (c) The particle density of cerium (blue) and oxygen (orange) atoms at a grain boundary in cerium oxide. (d) The electrostatic potential at a cerium oxide grain boundary. (e) The center of mass of cerium (blue) and oxygen (orange) atoms in two dimensions, at a grain boundary in cerium oxide.
</figcaption>
</p>

# `polypy`

`polypy` is a Python module for analysing molecular dynamics and Monte Carlo trajectories generated from the DL_POLY [@Smith] and DL_MONTE [@Purton] codes. The code reads DL_POLY HISTORY and CONFIG files, DL_MONTE ARCHIVE files and stores the data in a `polypy.read.Trajectory` object which is then used by the various data analysis modules.

The `polypy.density.Density` module generates a three dimensional grid and counts the number of times a given atom spends at each grid point during the simulation. This is then used to generate the particle density of a given atom in one and two dimensions. From here, the charge density in one and two dimensions, the electric field in one dimension and electrostatic potential in one dimension can be calculated using the `polypy.analysis` module. 

The `polypy.msd` module performs a mean squared displacement calculation. From the mean squared displacement, the three, two and one dimensional diffusion coefficient and ionic coefficient can be calculated. 

A module allowing easy generation of publication plots from the calculated data is available. The outputs are returned in a sensible form, allowing further manipulation and plotting.
The repository contains examples of the core functionality as well as tutorials, implemented in Jupyter notebooks to explain the full theory. Furthermore, a detailed description of theory is also available within the documentation. `polypy` is aimed towards theoretical solid state physicists and chemists who have a basic familiarity with Python, although the examples contained in the repository are designed to help less experienced users make use of the code.

# Acknowledgements
  
This package was written during a PhD project that was funded by AWE and EPSRC (EP/R010366/1). The `polypy` software package was developed to analyse data generated using the Balena HPC facility at the University of Bath and the ARCHER UK National Supercomputing Service (http://www.archer.ac.uk) via our membership of the UK's HEC Ma-terials Chemistry Consortium funded by EPSRC (EP/L000202).The author would like to thank Andrew R. McCluskey, Benjamin Morgan, Marco Molinari, James Grant and Stephen C. Parker for their help and guidance during this PhD project.

# References