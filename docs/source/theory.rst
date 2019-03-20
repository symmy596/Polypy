Theory
======

`polypy` is a Python module to analyse DLPOLY and DLMONTE trajectory files. Before using this code you will need to generate the relevant data. `polypy` is aimed at the solid state community and there are a wide range of applications. 

Charge Density
--------------

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


Ionic Conductivity
------------------


