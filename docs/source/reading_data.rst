Tutorial 1 - Reading data
-------------------------

The HISTORY, ARCHIVE and CONFIG classes expects two things, the filename
of the history file and a list of atoms to read. They will return a
``polypy.read.Trajectory`` object, which stores the the atom labels
(``Trajectory.atom_list``), datatype (``Trajectory.datatype``),
cartesian coordinates (``Trajectory.cartesian_coordinates``), fractiona
coordinates (``Trajectory.fractional_coordinates``), reciprocal lattice
vectors (``Trajectory.reciprocal_lv``), lattice vectors
(``Trajectory.lv``) cell lengths (``Trajectory.cell_lengths``), total
atoms in the file (``Trajectory.atoms_in_history``), timesteps
(``Trajectory.timesteps``), total atoms per timestep
(``Trajectory.total_atoms``).

HISTORY Files
~~~~~~~~~~~~~

.. code:: ipython3

    from polypy import read as rd

.. code:: ipython3

    history = rd.History("../example_data/HISTORY_CaF2", ["CA", "F"])

.. code:: ipython3

    print(history.trajectory.fractional_trajectory)


.. parsed-literal::

    [[0.5170937  0.51658126 0.51643485]
     [0.51658126 0.61669107 0.61654466]
     [0.61669107 0.51658126 0.61691069]
     ...
     [0.46866197 0.25395423 0.58485915]
     [0.37035211 0.58795775 0.45221831]
     [0.36552817 0.48637324 0.17484859]]


.. code:: ipython3

    print(history.trajectory.timesteps)


.. parsed-literal::

    500


.. code:: ipython3

    print(history.trajectory.atoms_in_history)


.. parsed-literal::

    750000


.. code:: ipython3

    print(history.trajectory.total_atoms)


.. parsed-literal::

    1500


It is often necessary to remove the equilibriation timesteps from the
simulation. This can be accomlished with the remove_initial_timesteps
method to remove timesteps at the start of the simulation and the
remove_final_timesteps, to remove timesteps at the end of the
simulation.

.. code:: ipython3

    new_history = history.trajectory.remove_initial_timesteps(10)
    print(new_history.timesteps)


.. parsed-literal::

    490


.. code:: ipython3

    new_history = new_history.remove_final_timesteps(10)
    print(new_history.timesteps)


.. parsed-literal::

    480


It is possible to return the trajectory for a single timestep within the
history file or to return the trajectory for a single atom.

.. code:: ipython3

    config_ca = history.trajectory.get_atom("CA")
    
    print(config_ca.fractional_trajectory)


.. parsed-literal::

    [[0.5170937  0.51658126 0.51643485]
     [0.51658126 0.61669107 0.61654466]
     [0.61669107 0.51658126 0.61691069]
     ...
     [0.31458099 0.41869718 0.41764085]
     [0.42742958 0.32461268 0.42507042]
     [0.42485915 0.42183099 0.31564789]]


.. code:: ipython3

    config_1 = history.trajectory.get_config(1)
    
    print(config_1.fractional_trajectory)


.. parsed-literal::

    [[0.53227339 0.51016082 0.50950292]
     [0.52116228 0.62894737 0.61761696]
     [0.62240497 0.50526316 0.6056652 ]
     ...
     [0.39444444 0.44974415 0.45102339]
     [0.45599415 0.37865497 0.39890351]
     [0.36343202 0.49309211 0.3690424 ]]


CONFIG Files
~~~~~~~~~~~~

.. code:: ipython3

    config = rd.Config("../example_data/CONFIG", ["CA", "F"])

.. code:: ipython3

    print(config.trajectory.fractional_trajectory)


.. parsed-literal::

    [[0.51666667 0.51666667 0.51666667]
     [0.51666667 0.61666667 0.61666667]
     [0.61666667 0.51666667 0.61666667]
     ...
     [0.36666667 0.46666667 0.46666667]
     [0.46666667 0.36666667 0.36666667]
     [0.36666667 0.46666667 0.36666667]]


DLMONTE
~~~~~~~

.. code:: ipython3

    archive = rd.Archive("../example_data/ARCHIVE_Short", ["AL"])

.. code:: ipython3

    print(archive.trajectory.timesteps)


.. parsed-literal::

    1000


