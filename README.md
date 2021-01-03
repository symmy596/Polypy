<p align="center"> 
<img src="https://github.com/symmy596/Polypy/blob/master/docs/source/Figures/polypy_1.png?raw=true"/>
</p>

[![PyPI version](https://badge.fury.io/py/polypy.svg)](https://badge.fury.io/py/polypy)
[![Build Status](https://travis-ci.com/symmy596/PolyPy.svg?branch=master)](https://travis-ci.com/symmy596/PolyPy)
[![Build status](https://ci.appveyor.com/api/projects/status/eo426m99lmkbh5rx?svg=true)](https://ci.appveyor.com/project/symmy596/polypy)
[![Documentation Status](https://readthedocs.org/projects/polypy/badge/?version=latest)](https://polypy.readthedocs.io/en/latest/?badge=latest)
<a href='https://coveralls.io/github/symmy596/PolyPy?branch=master'><img src='https://coveralls.io/repos/github/symmy596/PolyPy/badge.svg?branch=master' alt='Coverage Status' /></a> 

This is the documentation for the open-source Python project - `polypy`,
A library designed to facilitate the analysis of [DL_POLY](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx) and [DL_MONTE](https://www.ccp5.ac.uk/DL_MONTE) calculations.
polypy is built on existing Python packages that those in the solid state physics/chemistry community should already be familiar with.
It is hoped that this tool will bring some benfits to the solid state community and facilitate data analysis and the generation of publication ready plots (powered by Matplotlib.)

The main features include:

1. **Method to analyse the number denisty of a given species in one and two dimensions.**  

   - Generate a plot of the total number of species in bins perpendicular to a specified direction.  
   - Generate a plot of the total number of species in cuboids parallel to a specified direction.  

2. **Method calculate the charge density from the number density.**  

   - Convert number densities of all species in bins perpendicular to a specified direction into the charge density.  

3. **Calculate the electric field and electrostatic potential from the charge density.**  

   - Solves the Poisson Boltzmann equation to convert the charge density into the electric field and the electrostatic potential.

4. **Calculate the diffusion coefficient for a given species from a mean squared displacement.**

   - Carries out a mean squared displacement on an MD trajectory.
   - Calculates the diffusion coefficient.
   - Uses the density analysis and the diffusion coefficient to calculate the ionic conductivity. 
   

<p align="center"> 
<img src="https://github.com/symmy596/PolyPy/blob/master/docs/source/Figures/Show_off.png?raw=true" width="100%"/>
<figcaption>(a) Mean squared displacement for calcuim fluoride. (b) System volume of calcium fluoride during a molecular dynamics simulation. (c) Cerium and oxygen number density at a grain boundary. (d) Electrostatic potential across a grain boundary of cerium oxide. (e) Center of mass of cerium (blue) and oxygen (orange) in a cerium oxide grain boundary in two dimensions.</figcaption>
</p>


The code has been developed to analyse DL_POLY and DL_MONTE calculations however other codes can be incorporated if there is user demand. Other formats, such as pdb or xyz can be converted to `DL_POLY` format with codes such as [atomsk](https://atomsk.univ-lille.fr/) and then analysed with `polypy`. Users are welcome to increase the file coverage by adding a reading function for a different format. This can be accomplished by adding to the `read` module which has a class for each unique file type that converts it to a `polypy.read.trajectory` object. 

`polypy` was developed during a PhD project and as such the functionality focuses on the research questions encountered during that project, which we should clarify are wide ranging. Code contributions aimed at expanding the code to new of problems are encouraged.

`polypy` is free to use.

## Usage

A full list of examples can be found in the examples folder of the git repository, these include both the Python scripts and jupyter notebook tutorials which combine the full theory with code examples. It should be noted however that DL_POLY HISTORY files and DL_MONTE ARCHIVE files are sizable (1-5GB) and as such only short example trajectories are included in this repository. Notebooks are provided here to illustrate the theory but are not practicle.

## Installation

`polypy` is a Python 3 package and requires a typical scientific Python stack. Use of the tutorials requires Anaconda/Jupyter to be installed.

To build from source:

    pip install -r requirements.txt

    python setup.py build

    python setup.py install

Or alternatively install with pip

    pip install polypy

Using conda, 

    conda skeleton pypi polypy

    conda build polypy
    
    conda install --use-local polypy


### Tests

Tests can be run by typing:

    python setup.py test


in the root directory. 


### Documentation

To build the documentation from scratch
  
    cd docs
    make html

An online version of the documentation can be found [here](https://polypy.readthedocs.io/en/latest/index.html). The documentation contains an extensive explanation of the underlying theory, function documentation and tutorials. 


### License

`polypy` is made available under the MIT License.

### Detailed requirements

`polypy` is compatible with Python 3.5+ and relies on a number of open source Python packages, specifically:

- matplotlib
- numpy
- scipy
- coveralls
- coverage
- seaborn
- pandas
- jupyter
- nbsphinx
- jupyter-sphinx==0.2.4
- sphinx_rtd_theme

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
["fork and pull"](https://guides.github.com/activities/forking/)
workflow while core developers use branches in the main repository:

   1. First open an Issue to discuss the proposed contribution. This
      discussion might include how the changes fit polypy's scope and a
      general technical approach.
   2. Make your own project fork and implement the changes
      there. Please keep your code style compliant with PEP8.
   3. Open a pull request to merge the changes into the main
      project. A more detailed discussion can take place there before
      the changes are accepted.

For further information please contact Adam Symington, ars44@bath.ac.uk

## Future

Listed below are a series of useful additions that we would like to make to the codebase. Users are encouraged to fork the repository and work on any of these problems. Indeed, if functionality is not listed below you are more than welcome to add it. 

- RDF
- Diagonal slices
- Regional MSDs in a cube


## Author

* Adam R.Symington
  
## Research 

- [Defect segregation facilitates oxygen transport at fluorite UO2 grain boundaries](https://royalsocietypublishing.org/doi/full/10.1098/rsta.2019.0026)
- [The role of dopant segregation on the oxygen vacancy distribution and oxygen diffusion in CeO2 grain boundaries](https://iopscience.iop.org/article/10.1088/2515-7655/ab28b5/meta)
- [Quantifying the impact of disorder on Li-ion and Na-ion transport in perovskite titanate solid electrolytes for solid-state batteries](https://pubs.rsc.org/en/content/articlehtml/2020/ta/d0ta05343k)


## Acknowledgements
  
This package was written during a PhD project that was funded by AWE and EPSRC (EP/R010366/1). The `polypy` software package was developed to analyse data generated using the Balena HPC facility at the University of Bath. The author would like to thank Andrew R. McCluskey, Benjamin Morgan, Marco Molinari and Stephen C. Parker for their help and guidance during this PhD project.

