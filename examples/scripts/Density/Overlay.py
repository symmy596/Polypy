from polypy import read as rd
from polypy import density as Dens
from polypy import utils as ut
from polypy import write as wr
import numpy as np


data = rd.read_history("../example_data/HISTORY", ["CA", "F"])

ca_density = Dens.Density(data, atom_type="CA")
f_density = Dens.Density(data, atom_type="F")

cx, cy, cz, cy2 = ca_density.one_and_two_dimension_overlay(box=0.1)
fx, fy, fz, fy2 = f_density.one_and_two_dimension_overlay(box=0.1)

wr.combined_density_plot(cx, cy, cz, cy2)
wr.combined_density_plot(fx, fy, fz, fy2)