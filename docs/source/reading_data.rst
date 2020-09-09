Reading Data
============

The HISTORY, ARCHIVE and CONFIG classes expects two things, the filename of the history file and a list of atoms to read. They will return a `polypy.read.Trajectory` object, which stores the the atom labels (`Trajectory.atom_labels`), datatype (`Trajectory.data_type`), cartesian coordinates (`Trajectory.cartesian_coordinates`), fractiona coordinates (`Trajectory.fractional_coordinates`), reciprocal lattice vectors (`Trajectory.reciprocal_lv`), lattice vectors (`Trajectory.lv`) cell lengths (`Trajectory.cell_lengths`), total atoms in the file (`Trajectory.atoms_in_history`), timesteps (`Trajectory.timesteps`), total atoms per timestep (`Trajectory.total_atoms`).

HISTROY Files
~~~~~~~~~~~~~

.. code-block:: python

    from polypy import read as rd

    history = rd.History("../example_data/HISTORY", ["CA", "F"])    
    print(history.trajectory.fractional_trajectory)

It is possible to return the trajectory for a single timestep within the history file or to return the trajectory for a single atom.

.. code-block:: python

    config_ca = history.trajectory.get_atom("CA")

    print(config_ca.fractional_trajectory)  


.. code-block:: python

    config_1 = history.trajectory.get_config(1)

    print(config_1.fractional_trajectory)


CONFIG Files
~~~~~~~~~~~~

.. code-block:: python

    config = rd.Config("../example_data/CONFIG", ["CA", "F"])
    print(config.trajectory.fractional_trajectory)

ARCHIVE Files
~~~~~~~~~~~~~

.. code-block:: python

    archive = rd.Archive("../example_data/ARCHIVE", ["CE"])
    print(archive.trajectory.fractional_trajectory)
