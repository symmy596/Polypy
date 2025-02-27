Mean Squared Displacement MSD
=============================

Mean Squared Displacement (MSD)
-------------------------------

Molecules in liquds, gases and solids do not stay in the same place and
move constantly. Think about a drop of dye in a glass of water, as time
passes the dye distributes throughout the water. This process is called
diffusion and is common throughout nature and an incredibly relevant
property for materials scientists who work on things like batteries.

Using the dye as an example, the motion of a dye molecule is not simple.
As it moves it is jostled by collisions with other molecules, preventing
it from moving in a straight path. If the path is examined in close
detail, it will be seen to be a good approximation to a random walk. In
mathmatics a random walk is a series of steps, each taken in a random
direction. This was analysed by Albert Einstein in a study of Brownian
motion and he showed that the mean square of the distance travelled by a
particle following a random walk is proportional to the time elapsed.

.. math::
    \Big \langle r_{i}^{2} \big \rangle = 6 D_t + C 

where 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big>

where :math:`\Big \langle r^2 \big \rangle` is the mean squared
distance, t is time, :math:`D_t` is the diffusion rate and C is a
constant. If :math:`\Big \langle r_{i}^{2} \big \rangle` is plotted as a
function of time, the gradient of the curve obtained is equal to 6 times
the self-diffusion coefficient of particle i. The state of the matter
effects the shape of the MSD plot, solids, where little to no diffusion
is occuring, has a flat MSD profile. In a liquid however, the particles
diffusion randomly and the gradient of the curve is proportional to the
diffusion coefficient.

What is the mean squared displacement
-------------------------------------

Going back to the example of the dye in water, lets assume for the sake
of simplicity that we are in one dimension. Each step can either be
forwards or backwards and we cannot predict which. From a given starting
position, what distance is our dye molecule likely to travel after 1000
steps? This can be determined simply by adding together the steps,
taking into account the fact that steps backwards subtract from the
total, while steps forward add to the total. Since both forward and
backward steps are equally probable, we come to the surprising
conclusion that the probable distance travelled sums up to zero.

By adding the square of the distance we will always be adding positive
numbers to our total which now increases linearly with time. Based upon
equation 1 it should now be clear that a plot of
:math:`\Big \langle r_{i}^{2} \big \rangle` vs time with produce a line,
the gradient of which is equal to 6D. Giving us direct access to the
diffusion coefficient of the system.

.. code:: ipython3

    from polypy import read as rd
    from polypy.msd import MSD 
    from polypy.msd import RegionalMSD 
    from polypy import analysis
    from polypy import utils as ut
    from polypy import plotting
    import numpy as np
    import matplotlib.pyplot as plt

This example will use a short (50,000 steps), pre-prepared trajectory of
bulk :math:`CaF_2`. In reality we probably want a considerably longer
simulation (~10,000,000 steps). Such simulations generate huge files
(5GB) and the analysis would take too long for this tutorial.

The first step is to read the history file to generate the data. The
``HISTORY`` class expects two things, the filename of the history file
and a list of atoms to read. It will return a ``polypy.read.Trajectory``
object, which stores the the atom labels (``Trajectory.atom_labels``),
datatype (``Trajectory.data_type``), cartesian coordinates
(``Trajectory.cartesian_coordinates``), fractiona coordinates
(``Trajectory.fractional_coordinates``), reciprocal lattice vectors
(``Trajectory.reciprocal_lv``), lattice vectors (``Trajectory.lv``) cell
lengths (``Trajectory.cell_lengths``), total atoms in the file
(``Trajectory.atoms_in_history``), timesteps (``Trajectory.timesteps``),
total atoms per timestep (``Trajectory.total_atoms``).

.. code:: ipython3

    history_caf2 = rd.History("../example_data/HISTORY_CaF2", ["F"])

Once the data has been read into the code the MSD calculation can be
performed using the ``MSD`` class. The code will return a
``polypy.MSD.MSDContainer`` object, which contains the MSD information.

.. code:: ipython3

    f_msd = MSD(history_caf2.trajectory, sweeps=2)
    
    output = f_msd.msd()

.. code:: ipython3

    ax = plotting.msd_plot(output)
    
    plt.show()



.. image:: Figures/output_6_0.png

MSD calculations require a large number of statistics to be considered representative. A full msd will use every single frame of the trajectory as a starting point and effectively do a seperate msd from each starting point, these are then averaged to give the final result.  An MSD is technically an ensemble average over all sweeps and number of particles. 
The sweeps paramter is used to control the number of frames that are used as starting points in the calculation. For simulations with lots of diffusion events, a smaller number will be sufficient whereas simulations with a small number of diffusion events will require a larger number. 

.. code:: ipython3

    f_msd = MSD(history_caf2.trajectory, sweeps=10)
    
    output = f_msd.msd()
    
    ax = plotting.msd_plot(output)
    plt.show()



.. image:: Figures/output_7_0.png


.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.6078332646337548
    One Dimensional Diffusion Coefficient in X 1.6045620180115865
    One Dimensional Diffusion Coefficient in Y 1.6856414148385679
    One Dimensional Diffusion Coefficient in Z 1.5332963610511103


Note:
An MSD is supposed to be linear only after a ballistic regime and it usually lacks statistics for longer times. Thus the linear fit to extract the slope and thus the diffusion coefficient should be done on a portion of the MSD only.
This can be accomplished using the `exclude_initial` and `exclude_final` parameters

.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.5912662736409342
    One Dimensional Diffusion Coefficient in X 1.5862517497696607
    One Dimensional Diffusion Coefficient in Y 1.6753802400942055
    One Dimensional Diffusion Coefficient in Z 1.5121668310589353


Arrhenius
---------

It is then possible to take diffusion coefficients, calculated over a
large temperature range and, using the Arrhenius equation calculate the
activation energy for diffusion. Common sense and chemical intuition
suggest that the higher the temperature, the faster a given chemical
reaction will proceed. Quantitatively this relationship between the rate
a reaction proceeds and its temperature is determined by the Arrhenius
Equation. At higher temperatures, the probability that two molecules
will collide is higher. This higher collision rate results in a higher
kinetic energy, which has an effect on the activation energy of the
reaction. The activation energy is the amount of energy required to
ensure that a reaction happens.

.. math::
    k = A * e^{(-Ea / RT)}

where k is the rate coefficient, A is a constant, Ea is the activation
energy, R is the universal gas constant, and T is the temperature (in
kelvin).

Ionic Conductivity
------------------

Usefully, as we have the diffusion coefficient, the number of particles
(charge carriers) and the ability to calculate the volume, we can
convert this data into the ionic conductivity and then the resistance.

.. math::
    \sigma = \frac{D C_F e^2}{k_B T} 


where :math:`\sigma` is the ionic conductivity, D is the diffusion
coefficient, :math:`C_F` is the concentration of charge carriers, which
in this case if F ions, :math:`e^2` is the charge of the diffusing
species, :math:`k_B` is the Boltzmann constant and T is the temperature.

The resitance can then be calculated according to

.. math::
    \Omega = \frac{1}{\sigma} 

So the first step is to calculate the volume, the system voume module
will do this from the given data.

.. code:: ipython3

    volume, step = analysis.system_volume(history_caf2.trajectory)
    average_volume = np.mean(volume[:50])

The number of charge carriers is just the total number of atoms.

.. code:: ipython3

    sigma = analysis.conductivity(history_caf2.trajectory.total_atoms, 
                            average_volume, 
                            output.xyz_diffusion_coefficient(), 
                            1500, 1)

.. code:: ipython3

    print("Ionic Conductivity :", sigma)


.. parsed-literal::

    Ionic Conductivity : 0.0008752727736501591


.. code:: ipython3

    print("Resistivity :", (1 / sigma)) 


.. parsed-literal::

    Resistivity : 1142.5009781004494


Simulation Length
-----------------

It is important to consider the lenght of your simulation (Number of
steps). The above examples use a short trajectory but it is at a
sufficient temperature that there are enough diffusion events to get a
good MSD plot. The following example is of a very short simulation, you
will hopefully note that the MSD plot is clearly not converged.

.. code:: ipython3

    history_short = rd.History("../example_data/HISTORY_short", atom_list=["F"])

.. code:: ipython3

    f_msd_short = MSD(history_short.trajectory, sweeps=2)
    
    output = f_msd_short.msd()

.. code:: ipython3

    ax = plotting.msd_plot(output)
    plt.show()



.. image:: Figures/output_19_0.png


.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.58656319093229
    One Dimensional Diffusion Coefficient in X 1.5739020833099904
    One Dimensional Diffusion Coefficient in Y 1.630216356788139
    One Dimensional Diffusion Coefficient in Z 1.5555711326987387


Amusingly, this actually does not seem to have a huge effect on the
diffusion coefficient compared to the longer simulation. However these
trajectories are from a CaF\ :math:`_2` simulation at 1500 K and there
are thus a large number of diffusion events in the short time frame.

State of Matter
---------------

It is possible to identify the phase of matter from the MSD plot.

.. raw:: html

   <center>

 Figure 1. The anticipated MSD form for each state of matter.

.. raw:: html

   </center>

The fluorine diffusion discussed already clearly shows that the fluorine
sub lattice has melted and the diffusion is liquid like. Whereas,
carrying out the same analysis on the calcium sub lattice shows that
while the fluorine sub lattice has melted, the Calcium sub lattice is
still behaving like a solid.

.. code:: ipython3

    f_msd = MSD(history_caf2.trajectory, sweeps=2)
    
    output = f_msd.msd()
    
    ax = plotting.msd_plot(output)
    plt.show()



.. image:: Figures/output_23_0.png


Regional MSD Calculations
-------------------------

Often in solid state chemistry simulations involve defects, both
structural e.g. grain boundaries, dislocations and surface, and chemical
e.g. point defects. It is important to try and isolate the contributions
of these defects to the overall properties. Regarding diffusion, it
could be imagined that a certain region within a structure will have
different properties compared with the stoichiometric bulk, e.g. a grain
boundary vs the grains, or the surface vs the bulk. ``polypy`` has the
capability to isolate trajectories that pass within certain regions of a
structure and thus calculate a diffusion coefficient for those regions.

In this example we will calculate the diffusion coefficient in a box
between -5.0 and 5.0 in the dimension of the first lattice vector.

.. code:: ipython3

    f_msd = RegionalMSD(history_caf2.trajectory, -5, 5, dimension="x")
    output = f_msd.analyse_trajectory()

.. code:: ipython3

    ax = plotting.msd_plot(output)
    
    plt.show()



.. image:: Figures/output_26.png


.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.597047044241002
    One Dimensional Diffusion Coefficient in X 1.6120172452124801
    One Dimensional Diffusion Coefficient in Y 1.671268658071343
    One Dimensional Diffusion Coefficient in Z 1.5078552294391845


DLMONTE
^^^^^^^

.. code:: ipython3

    archive = rd.Archive("../example_data/ARCHIVE_LLZO", atom_list=["O"])

.. code:: ipython3

    f_msd = MSD(archive.trajectory, sweeps=2)


::


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-20-2e636209fda5> in <module>
    ----> 1 f_msd = MSD(archive.trajectory, sweeps=2)
    

    /opt/anaconda3/lib/python3.7/site-packages/polypy-0.7-py3.7.egg/polypy/msd.py in __init__(self, data, sweeps)
        153             raise ValueError("ERROR: MSD can only handle one atom type. Exiting")
        154         if data.data_type == "DL_MONTE ARCHIVE":
    --> 155             raise ValueError("DLMONTE simulations are not time resolved")
        156         self.distances = []
        157         self.msd_information = MSDContainer()


    ValueError: DLMONTE simulations are not time resolved

