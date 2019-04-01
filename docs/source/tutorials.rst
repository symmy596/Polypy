Tutorials
=========

These tutorials are replicated in jupyter notebook form and contained within examples. The accompanying python scripts are also included here. All of code examples within these tutorials can be found in `examples/Scripts <https://github.com/symmy596/PolyPy/tree/master/examples/scripts>`_.

All tutorials use fluorite :math:`CaF_2` as an example. 

Reading Data
------------

Density Analysis
----------------

Understanding the particle density of a material is incredibly useful when studying things like defect segregation and atomic structure. Consider a system with a grain boundary, it may be interesting to know the structure changes at the boundary, e.g is there an increase or decrease in the amount of a certain species at the grain boundary and does this inform you about aby segregation behaviour? Another example involves a surface with a layer of water

The density module of polypy allows the particle density to be evaluated in one and two dimensions, this can then be converted into a charge density and (in one dimension) the electric field and electrostatic potential.

.. toctree::
    :maxdepth: 1

    density_tutorial_1.rst
    density_tutorial_2.rst
    density_tutorial_3.rst

.. toctree::
    :maxdepth: 1

    msd_tutorial.rst
    
Mean Squared Displacement (MSD)
-------------------------------

.. toctree::
    :maxdepth: 1

    msd_tutorial.rst
    