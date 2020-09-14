Reading Data
============

The HISTORY, ARCHIVE and CONFIG classes expects two things, the filename of the history file and a list of atoms to read. They will return a `polypy.read.Trajectory` object, which stores;

- The the atom labels (`Trajectory.atom_labels`)
- The datatype (`Trajectory.data_type`)
- The cartesian coordinates (`Trajectory.cartesian_coordinates`)
- The fractional coordinates (`Trajectory.fractional_coordinates`)
- The reciprocal lattice vectors (`Trajectory.reciprocal_lv`)
- The lattice vectors (`Trajectory.lv`)
- The cell lengths (`Trajectory.cell_lengths`)
- The total atoms in the file (`Trajectory.atoms_in_history`)
- The timesteps (`Trajectory.timesteps`)
- The total atoms per timestep (`Trajectory.total_atoms`)

HISTROY Files
~~~~~~~~~~~~~

.. code-block:: python

    from polypy import read as rd

    history = rd.History("../example_data/HISTORY", ["CA", "F"])    

It is possible to return the trajectory for a single timestep within the history file or to return the trajectory for a single atom.

.. code-block:: python

    config_ca = history.trajectory.get_atom("CA")


.. code-block:: python

    config_1 = history.trajectory.get_config(1)


CONFIG Files
~~~~~~~~~~~~

.. code-block:: python

    config = rd.Config("../example_data/CONFIG", ["CA", "F"])

ARCHIVE Files
~~~~~~~~~~~~~

.. code-block:: python

    archive = rd.Archive("../example_data/ARCHIVE", ["CE"])
