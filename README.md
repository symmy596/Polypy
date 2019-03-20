# PolyPy

[![PyPI version](https://badge.fury.io/py/polypy.svg)](https://badge.fury.io/py/polypy) [![Build Status](https://travis-ci.com/symmy596/PolyPy.svg?branch=master)](https://travis-ci.com/symmy596/PolyPy) [![Coverage Status](https://coveralls.io/repos/github/symmy596/PolyPy/badge.svg?branch=master)](https://coveralls.io/github/symmy596/PolyPy?branch=master) [![Documentation Status](https://readthedocs.org/projects/polypy/badge/?version=latest)](https://polypy.readthedocs.io/en/latest/?badge=latest)


This is the documentation for the open-source Python project - `polypy`.
A library designed to facilitate the analysis of [DL_POLY](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx) and [DL_MONTE](https://www.ccp5.ac.uk/DL_MONTE) calculations.
polypy is built on existing Python packages that those in the solid state physics/chemistry community should already be familiar with.
It is hoped that this tool will bring some benfits to the solid state community and facilitate data analysis and the generation of publication ready plots (powered by Matplotlib.)

The main features include:

1. **Method to analyse the number denisty of a given species in one and two dimensions.**  

   - Generate a plot of the total number of species in bins perpendicular to a specified direction.  
   - Generate a plot of the total number of species in cuboids parallel to a specified direction.  
   - Detertmine the total number of species within a specific area of the system.

2. **Method calculate the charge density from the number density.**  

   - Convert number densities of all species in bins perpendicular to a specified direction into the charge density.  

3. **Calculate the electric field and electrostatic potential from the charge density.**  

   - Solves the Poisson Boltzmann equation to convert the charge density into the electric field and the electrostatic potential.

4. **Calculate the diffusion coefficient for a given species from a mean squared displacement.**

   - Carries out a mean squared displacement on an MD trajectory.
   - Calculates the diffusion coefficient.
   - Uses the density analysis and the diffusion coefficient to calculate the ionic conductivity. 

The code has been developed to analyse DL_POLY and DL_MONTE calculations however other codes can be incorporated if there is user demand.
`polypy` was developed during a PhD project and as such the functionality focuses on the research questions encountered during that project, which we should clarify
are wide ranging. Code contributions aimed at expanding the code to new of problems are encouraged.

`polypy` is free to use.

## Usage

A full list of examples can be found in the examples folder of the git repository, these include both the Python scripts and jupyter notebook tutorials which combine the full theory with code examples. It should be noted however that DL_POLY HISTORY files and DL_MONTE ARCHIVE files are sizable (1-5GB) and as such we would discourage the use of notebooks and encourage using `polypy` on the HPC resource used to generate the data. Notebooks are provided here to illustrate the theory but are not practicle.

## Installation

`polypy` is a Python 3 package and requires a typical scientific Python stack. Use of the tutorials requires Anaconda/Jupyter to be installed.

To build from source:

    pip install -r requirements.txt

    python setup.py build

    python setup.py install

    python setup.py test

Or alternatively install with pip

    pip install polypy


### Documentation

To build the documentation from scratch
  
    cd docs
    make html

### License

`polypy` is made available under the MIT License.

### Detailed requirements

`polypy` is compatible with Python 3.5+ and relies on a number of open source Python packages, specifically:

- Numpy
- Scipy
- Matplotlib

## Contributing

### Contact

If you have questions regarding any aspect of the software then please get in touch with the developer Adam Symington via email - ars44@bath.ac.uk.
Alternatively you can create an issue on the [Issue Tracker](https://github.com/symmy596/PolyPy/issues).

### Bugs

There may be bugs. If you think you've caught one, please report it on the [Issue Tracker](https://github.com/symmy596/PolyPy/issues).
This is also the place to propose new ideas for features or ask questions about the design of `polypy`. Poor documentation is considered a bug
so feel free to request improvements.

### Code contributions

We welcome help in improving and extending the package. This is managed through Github pull requests; for external contributions we prefer the
["fork and pull"](https://guides.github.com/activities/forking/)__
workflow while core developers use branches in the main repository:

   1. First open an Issue to discuss the proposed contribution. This
      discussion might include how the changes fit surfinpy's scope and a
      general technical approach.
   2. Make your own project fork and implement the changes
      there. Please keep your code style compliant with PEP8.
   3. Open a pull request to merge the changes into the main
      project. A more detailed discussion can take place there before
      the changes are accepted.

For further information please contact Adam Symington, ars44@bath.ac.uk

## Author

* Adam R.Symington
  
## Acknowledgements
  
* [Prof Stephen C.Parker](http://people.bath.ac.uk/chsscp/) - (Bath University)
