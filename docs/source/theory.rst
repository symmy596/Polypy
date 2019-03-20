Theory
======

`polypy` is a Python module to analyse DLPOLY and DLMONTE trajectory files. Before using this code you will need to generate the relevant data. `polypy` is aimed at the solid state community and there are a wide range of applications. 

Charge Density
--------------

Using the density module you can calculate the number density of atoms

.. math::
    \rho_{q}(z) = \sum_{i} q_{i} \rho_{i}(z)

where :math:`\rho_{i}` is the density of atom i and :math:`q_{i}` is its charge.    

Electric Field
--------------

.. math::
    E(z) = \frac{1}{- \epsilon_{0}} \int_{z_{0}}^{z} \rho_{q}(z')dz'

Electrostatic Potential
-----------------------

.. math::
    \Delta_{\psi}(z) = \int_{z_{0}}^{z} E(z')dz'



Mean Squared Displacement
-------------------------

As an atom travels it is interrupted by collisions with other atoms and this prevents it from travelling in a straight line. In the end the particle moves in a way resembling a random walk. In mathmatics a random walk is a series of steps, where each step is taken in a random direction. Albert Einstein showed that the mean square of the distance travelled by an atom following a random walk is proportional to the time elapsed. This relationship can be written as 

.. math::
    <r_{i}^{2}(t)> = 2Dt

where 

.. math::
    <r_{i}^{2}> = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big> 


Therefore if :math:`\Big< | r_{i}(t) - r_{i}(0) |^2 \Big>` is plotted as a function of time, the gradient of the curve obtained is equal to 6 times the self-diffusion coefficient of particle i. 
The state of the matter effects the shape of the MSD plot, solids, where little to no diffusion is occuring, has a flat MSD profile. In a liquid however, the particles diffusion randomly and the gradient of the curve is proportional to the diffusion coefficient. 

.. image:: Figures/MSD_Theory.png
    :height: 300px
    :align: center


Ionic Conductivity
------------------


