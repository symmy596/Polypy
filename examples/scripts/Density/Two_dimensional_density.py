from polypy import read as rd
from polypy import density as Dens
from polypy import utils as ut
from polypy import write as wr
import numpy as np


data = rd.read_history("../example_data/HISTORY", ["CA", "F"])

ca_density = Dens.Density(data, atom_type="CA")
f_density = Dens.Density(data, atom_type="F")

cx, cy, cz = ca_density.two_dimensional_density(box=0.1, direction="x")
fx, fy, fz = f_density.two_dimensional_density(box=0.1, direction="x")

wr.two_dimensional_density_plot(cx, cy, cz)
wr.two_dimensional_density_plot(fx, fy, fz)