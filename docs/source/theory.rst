Theory
======

`polypy` is a Python module to analyse DLPOLY and DLMONTE trajectory files. Before using this code you will need to generate the relevant data. `polypy` is aimed at the solid state community and there are a wide range of applications. 

Charge Density
--------------

Using the density module you can calculate the number density of atoms

.. math::
    \rho_{q}(z) = \sum_{i} q_{i} \rho_{i}(z)

where :math:`\rho_{i}` is the density of atom i and :math:`q_{i}` is its charge.    

Electric Field and Electrostatic Potential
------------------------------------------

The charge density can be converted into the electric field and the electrostatic potential.

The electric field is calculated according to 

.. math::
    E(z) = \frac{1}{- \epsilon_{0}} \int_{z_{0}}^{z} \rho_{q}(z')dz'

where :math:`\epsilon_{0}` is the permittivity of the vacuum and :math:`\rho_{q}` is the charge density.  
The electrostatic potential is calculated according to

.. math::
    \Delta_{\psi}(z) = \int_{z_{0}}^{z} E(z')dz'

Mean Squared Displacement
-------------------------

Molecules in liquds, gases and solids do not stay in the same place and move constantly. Think about a drop of dye in a glass of water, as time passes the dye distributes throughout the water. This process is called diffusion and is common throughout nature and an incredibly relevant property for materials scientists who work on things like batteries.  

Using the dye as an example, the motion of a dye molecule is not simple. As it moves it is jostled by collisions with other molecules, preventing it from moving in a straight path. If the path is examined in close detail, it will be seen to be a good approximation to a random walk. In mathmatics a random walk is a series of steps, each taken in a random direction. This was analysed by Albert Einstein in a study of Brownian motion and he showed that the mean square of the distance travelled by a particle following a random walk is proportional to the time elapsed. 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = 6 D_t + C 

where 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big>,


where :math:`\Big \langle r^2 \big \rangle` is the mean squared distance, t is time, :math:`D_t` is the diffusion rate and C is a constant. If :math:`\Big \langle r_{i}^{2} \big \rangle` is plotted as a function of time, the gradient of the curve obtained is equal to 6 times the self-diffusion coefficient of particle i. 

What is the mean squared displacement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Going back to the example of the dye in water, lets assume for the sake of simplicity that we are in one dimension. Each step can either be forwards or backwards and we cannot predict which. From a given starting position, what distance is our dye molecule likely to travel after 1000 steps? This can be determined simply by adding together the steps, taking into account the fact that steps backwards subtract from the total, while steps forward add to the total. Since both forward and backward steps are equally probable, we come to the surprising conclusion that the probable distance travelled sums up to zero.

By adding the square of the distance we will always be adding positive numbers to our total which now increases linearly with time. Based upon equation 1 it should now be clear that a plot of :math:`\Big \langle r_{i}^{2} \big \rangle` vs time with produce a line, the gradient of which is equal to 6D. Giving us direct access to the diffusion coefficient of the system. The state of the matter effects the shape of the MSD plot, solids, where little to no diffusion is occuring, has a flat MSD profile. In a liquid however, the particles diffusion randomly and the gradient of the curve is proportional to the diffusion coefficient. 

.. image:: Figures/States_of_Matter.png
    :height: 300px
    :align: center

The following example is for fluorine diffusion in :math:`CaF_2`.

.. image:: Figures/MSD_Theory.png
    :height: 300px
    :align: center


Ionic Conductivity
------------------

Usefully, as we have the diffusion coefficient, the number of particles (charge carriers) and the ability to calculate the volume, we can convert this data into the ionic conductivity and then the resistance. 

.. math::
    \sigma = \frac{D C_F e^2}{k_B T} 

where :math:`\sigma` is the ionic conductivity, D is the diffusion coefficient,:math:`C_F` is the concentration of charge carriers, which in this case if F ions, :math:`e^2` is the charge of the diffusing species, :math:`k_B` is the Boltzmann constant and T is the temperature. 

The resitance can then be calculated according to 

.. math::
    \Omega = \frac{1}{\sigma} 


Arrhenius
---------

It is possible to calculate the diffusion coefficients over a large temperature range and then use the Arrhenius equation to calculate the activation energy for diffusion. Common sense and chemical intuition suggest that the higher the temperature, the faster a given chemical reaction will proceed. Quantitatively this relationship between the rate a reaction proceeds and its temperature is determined by the Arrhenius Equation. At higher temperatures, the probability that two molecules will collide is higher. This higher collision rate results in a higher kinetic energy, which has an effect on the activation energy of the reaction. The activation energy is the amount of energy required to ensure that a reaction happens.  
  
.. math::
    k = A e^{(-Ea / RT)}
  
where k is the rate coefficient, A is a constant, Ea is the activation energy, R is the universal gas constant, and T is the temperature (in kelvin).