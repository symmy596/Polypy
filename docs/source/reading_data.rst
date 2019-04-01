Reading Data
============

The first thing to do with polypy is to read the data, whether it be in the DL_POLY HISTORY or CONFIG format. Polypy will read these formats and return a dictionary

.. code-block:: python

    data = {'label': Unique atom names, 
            'trajectories': Atomic trajectories, 
            'lv':Lattice vectors, 
            'timesteps':Total number of timesteps, 
            'natoms': Total number of atoms}

HISTROY Files
~~~~~~~~~~~~~

.. code-block:: python

    from polypy import read as rd
    import numpy as np

    history = rd.read_history("../example_data/HISTORY", ["CA", "F"])
    
    print(history['label'])

| ['CA' 'CA' 'CA' ... 'F' 'F' 'F']

polypy returns an (timesteps * number of atoms) x 3 array containing the atomic coordinates of all selected atoms at all timesteps. It is possible to isolate the CONFIG at a specified timestep using the get_config function.

.. code-block:: python

    print(history['trajectories'].size)
    print(history['trajectories'].shape)

    config_1 = rd.get_config(history, timestep=0)

    print(config_1) 

| 2250000
| (750000, 3)
| [[-13.193 -13.207 -13.211]
| [-13.207 -10.472 -10.476]
| [-10.472 -13.207 -10.466]
| ...
| [ 10.014  12.748  12.745]
| [ 12.749  10.025  10.014]
| [ 10.009  12.756  10.007]]

The total number of timesteps are returned as well as the total number of atoms.

.. code-block:: python

    print(history['timesteps'])
    print(history['natoms'])

| 500
| 1500

For things like mean squared displacements and particle densities it is usually best to split the trajectories for each species up into individual objects e.g. An MSD for just the F atoms. This is done by the get_atom function. 

.. code-block:: python

    f_data = rd.get_atom(history, "F")

    print(f_data['timesteps'])
    print(f_data['natoms'])

| 500
| 1000

Finally, it is possible to isolate the coordinates of a specific timestep in the HISTORY trajectory. 

.. code-block:: python

    first_config = rd.get_config(history, 1)
    print(first_config)

| [[-12.797  -13.402  -13.42  ]
|  [-13.101  -10.152  -10.462 ]
|  [-10.331  -13.536  -10.789 ]
|  ...
|  [ 10.792   12.305   12.34  ]
|  [ 12.476   10.36    10.914 ]
|  [  9.9435  13.491   10.097 ]]


.. code-block:: python

    first_config = rd.get_config(history, 1)

    print(first_config)

| [[-12.797  -13.402  -13.42  ]
|  [-13.101  -10.152  -10.462 ]
|  [-10.331  -13.536  -10.789 ]
|  ...
|  [ 10.792   12.305   12.34  ]
|  [ 12.476   10.36    10.914 ]
|  [  9.9435  13.491   10.097 ]]

CONFIG Files
~~~~~~~~~~~~

It is also possible to read CONFIG files and split the data into individual species as with the history file. 

.. code-block:: python

    config = rd.read_config("../example_data/CONFIG", ["CA", "F"])

    print(config['timesteps'])
    print(config['natoms'])

| 1
| 1500

.. code-block:: python

    f_data = rd.get_atom(config, "F")
    print(f_data['timesteps'])
    print(f_data['natoms'])

| 1
| 1000

ARCHIVE Files
~~~~~~~~~~~~~

.. code-block:: python

    archive = rd.read_archive("../example_data/ARCHIVE", ["CE", "O"])

    print(archive['timesteps'])
    print(archive['natoms'])

| 100
| 8640