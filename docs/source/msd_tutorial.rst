Mean Squared Displacement MSD
=============================

Molecules in liquds, gases and solids do not stay in the same place and move constantly. Think about a drop of dye in a glass of water, as time passes the dye distributes throughout the water. This process is called diffusion and is common throughout nature and an incredibly relevant property for materials scientists who work on things like batteries.  

Using the dye as an example, the motion of a dye molecule is not simple. As it moves it is jostled by collisions with other molecules, preventing it from moving in a straight path. If the path is examined in close detail, it will be seen to be a good approximation to a random walk. In mathmatics a random walk is a series of steps, each taken in a random direction. This was analysed by Albert Einstein in a study of Brownian motion and he showed that the mean square of the distance travelled by a particle following a random walk is proportional to the time elapsed. 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = 6 D_t + C 

where 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big>

where :math:`\Big \langle r^2 \big \rangle` is the mean squared distance, t is time, :math:`D_t` is the diffusion rate and C is a constant. If :math:`\Big \langle r_{i}^{2} \big \rangle` is plotted as a function of time, the gradient of the curve obtained is equal to 6 times the self-diffusion coefficient of particle i. 
The state of the matter effects the shape of the MSD plot, solids, where little to no diffusion is occuring, has a flat MSD profile. In a liquid however, the particles diffusion randomly and the gradient of the curve is proportional to the diffusion coefficient. 

What is the mean squared displacement?
--------------------------------------

Going back to the example of the dye in water, lets assume for the sake of simplicity that we are in one dimension. Each step can either be forwards or backwards and we cannot predict which. From a given starting position, what distance is our dye molecule likely to travel after 1000 steps? This can be determined simply by adding together the steps, taking into account the fact that steps backwards subtract from the total, while steps forward add to the total. Since both forward and backward steps are equally probable, we come to the surprising conclusion that the probable distance travelled sums up to zero.

By adding the square of the distance we will always be adding positive numbers to our total which now increases linearly with time. Based upon equation 1 it should now be clear that a plot of :math:`\Big \langle r_{i}^{2} \big \rangle` vs time with produce a line, the gradient of which is equal to 6D. Giving us direct access to the diffusion coefficient of the system. 

Usage
~~~~~

.. code-block:: python

    from polypy import read as rd
    from polypy.msd import MSD 
    from polypy.msd import RegionalMSD 
    from polypy import analysis
    from polypy import utils as ut
    from polypy import plotting
    import numpy as np
    import matplotlib.pyplot as plt

This example will use a short (50,000 steps), pre-prepared trajectory of bulk $CaF_2$. In reality we probably want a considerably longer simulation (~10,000,000 steps). Such simulations generate huge files (5GB) and the analysis would take too long for this tutorial. 

The first step is to read the history file to generate the data. The `HISTORY` class expects two things, the filename of the history file and a list of atoms to read. It will return a `polypy.read.Trajectory` object, which stores the the atom labels (`Trajectory.atom_labels`), datatype (`Trajectory.data_type`), cartesian coordinates (`Trajectory.cartesian_coordinates`), fractiona coordinates (`Trajectory.fractional_coordinates`), reciprocal lattice vectors (`Trajectory.reciprocal_lv`), lattice vectors (`Trajectory.lv`) cell lengths (`Trajectory.cell_lengths`), total atoms in the file (`Trajectory.atoms_in_history`), timesteps (`Trajectory.timesteps`), total atoms per timestep (`Trajectory.total_atoms`).

.. code-block:: python

    history = rd.History("../example_data/HISTORY", ["F"])

Once the data has been read into the code the MSD calculation can be performed.

.. code-block:: python

    f_msd = MSD(history.trajectory, sweeps=2)

    output = f_msd.msd()

    ax = plotting.msd_plot(output)
    plt.show()

.. image:: Figures/MSD_1.png
    :align: center

.. code-block:: python

    f_msd = MSD(history.trajectory, sweeps=10)

    output = f_msd.msd()

    ax = plotting.msd_plot(output)

    plt.show()

.. image:: Figures/MSD_2.png
    :align: center

Using the data the diffusion coefficient can then be calculated from the slopes. 

.. code-block:: python

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())

| Three Dimensional Diffusion Coefficient 1.6078332646337548
| One Dimensional Diffusion Coefficient in X 1.6045620180115865
| One Dimensional Diffusion Coefficient in Y 1.6856414148385679
| One Dimensional Diffusion Coefficient in Z 1.5332963610511103


Arrhenius
~~~~~~~~~

It is then possible to take diffusion coefficients, calculated over a large temperature range and, using the Arrhenius equation calculate the activation energy for diffusion. Common sense and chemical intuition suggest that the higher the temperature, the faster a given chemical reaction will proceed. Quantitatively this relationship between the rate a reaction proceeds and its temperature is determined by the Arrhenius Equation. At higher temperatures, the probability that two molecules will collide is higher. This higher collision rate results in a higher kinetic energy, which has an effect on the activation energy of the reaction. The activation energy is the amount of energy required to ensure that a reaction happens.  
  
.. math::
    k = A * e^{(-Ea / RT)}
  
where k is the rate coefficient, A is a constant, Ea is the activation energy, R is the universal gas constant, and T is the temperature (in kelvin).


Ionic Conductivity
~~~~~~~~~~~~~~~~~~

Usefully, as we have the diffusion coefficient, the number of particles (charge carriers) and the ability to calculate the volume, we can convert this data into the ionic conductivity and then the resistance. 

.. math::
    \sigma = \frac{D C_F e^2}{k_B T} 

where :math:`\sigma` is the ionic conductivity, D is the diffusion coefficient, :math:`C_F` is the concentration of charge carriers, which in this case if F ions, :math:`e^2` is the charge of the diffusing species, :math:`k_B` is the Boltzmann constant and T is the temperature. 

The resitance can then be calculated according to 

.. math::
    \Omega = \frac{1}{\sigma} 

So the first step is to calculate the volume, the system volume module will do this from the given data. 

.. code-block:: python

    volume, step = analysis.system_volume(history.trajectory)
    average_volume = np.mean(volume[:50])

The number of charge carriers is just the total number of atoms.

.. code-block:: python

    sigma = analysis.conductivity(history.trajectory.total_atoms, 
                                  average_volume, 
                                  output.xyz_diffusion_coefficient(), 
                                  1500)
    print("Ionic Conductivity :", sigma)

Ionic Conductivity : 0.0008752727736501591

.. code-block:: python

    print("Resistivity :", (1 / sigma))     

Resistivity : 1142.5009781004494

Simulation Length
~~~~~~~~~~~~~~~~~

It is important to consider the lenght of your simulation (Number of steps). The above examples use a short trajectory but it is at a sufficient temperature that there are enough diffusion events to get a good MSD plot. The following example is of a very short simulation, you will hopefully note that the MSD plot is clearly not converged.

.. code-block:: python

    data_short = rd.History("../example_data/HISTORY_short", atom_list=["F"])
    f_msd_short = MSD(data_short.trajectory, sweeps=2)

    output = f_msd_short.msd()

    ax = plotting.msd_plot(output)
    plt.show()

.. image:: Figures/MSD_3.png
    :align: center

.. code-block:: python

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())

| Three Dimensional Diffusion Coefficient 1.58656319093229
| One Dimensional Diffusion Coefficient in X 1.5739020833099904
| One Dimensional Diffusion Coefficient in Y 1.630216356788139
| One Dimensional Diffusion Coefficient in Z 1.5555711326987387


State of Matter
~~~~~~~~~~~~~~~

It is possible to identify the phase of matter from the MSD plot.

.. image:: Figures/States_of_Matter.png
    :height: 300px
    :align: center

The Fluorine diffusion discussed already clearly shows that the fluorine sub lattice has melted and the diffusion is liquid like. Whereas, carrying out the same analysis on the Calcium sub lattice shows that while the fluorine sub lattice has melted, the Calcium sub lattice is still behaving like a solid. 

.. code-block:: python

    history = rd.History("../example_data/HISTORY", ["CA"])

    f_msd = MSD(history.trajectory, sweeps=2)

    output = f_msd.msd()

    ax = plotting.msd_plot(output)
    plt.show()

.. image:: Figures/MSD_4.png
    :align: center

Regional MSD
~~~~~~~~~~~~

Often in solid state chemistry simulations involve defects, both structural e.g. grain boundaries, dislocations and surface, and chemical e.g. point defects. It is important to try and isolate the contributions of these defects to the overall properties. Regarding diffusion, it could be imagined that a certain region within a structure will have different properties compared with the stoichiometric bulk, e.g. a grain boundary vs the grains, or the surface vs the bulk. `polypy` has the capability to isolate trajectories that pass within certain regions of a structure and thus calculate a diffusion coefficient for those regions. 

.. code-block:: python

    history = rd.History("../example_data/HISTORY", atom_list=["F"])

    f_msd = RegionalMSD(history.trajectory, -5, 5)
    output = f_msd.analyse_trajectory()

    ax = plotting.msd_plot(output)

    plt.show()

.. image:: Figures/MSD_5.png
    :align: center


.. code-block:: python

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())

| Three Dimensional Diffusion Coefficient 1.597047044241002
| One Dimensional Diffusion Coefficient in X 1.6120172452124801
| One Dimensional Diffusion Coefficient in Y 1.671268658071343
| One Dimensional Diffusion Coefficient in Z 1.5078552294391845
