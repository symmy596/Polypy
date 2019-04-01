Density Within a Specified Region
=================================

The first feature of the density class allows the calculation of the particle density in one dimension of a given species. The total number of atoms between an upper and lower coordinate value e.g. the top and bottom coordinates of a slab, or edges of a grain boundary. 

The first step is to read the data. We want the data for both species so need to provide a list of the species.

.. code-block:: python

    from polypy import read as rd
    from polypy import density as Dens
    from polypy import utils as ut
    from polypy import write as wr
    import numpy as np  

    ["CA", "F"]

    data = rd.read_history("../example_data/HISTORY", ["CA", "F"])

The next step is to create the density object for both species.

.. code-block:: python

    ca_density = Dens.Density(data, atom_type="CA")
    f_density = Dens.Density(data, atom_type="F")


.. code-block:: python

    plane = ca_density.one_dimensional_density_sb(ul=5.0, ll=-5.0)
    print("Total Number of Ca Between -5.0 - 5.0   :", plane, " across ", data['timesteps'], "timesteps")
    print("Average Number of Ca Between -5.0 - 5.0 :", plane / data['timesteps'])

    plane = f_density.one_dimensional_density_sb(ul=5.0, ll=-5.0)
    print("Total Number of F Between -5.0 - 5.0   :", plane, " across ", data['timesteps'], "timesteps")
    print("Average Number of F Between -5.0 - 5.0 :", plane / data['timesteps'])

| Total Number of Ca Between -5.0 - 5.0   : 77999  across  500 timesteps
| Average Number of Ca Between -5.0 - 5.0 : 155.998
| Total Number of F Between -5.0 - 5.0   : 190004  across  500 timesteps
| Average Number of F Between -5.0 - 5.0 : 380.008