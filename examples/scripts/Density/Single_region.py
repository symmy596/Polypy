from polypy import read as rd
from polypy import density as Dens
from polypy import utils as ut
from polypy import write as wr
import numpy as np


data = rd.read_history("../example_data/HISTORY", ["CA", "F"])

ca_density = Dens.Density(data, atom_type="CA")
f_density = Dens.Density(data, atom_type="F")

plane = ca_density.one_dimensional_density_sb(ul=5.0, ll=-5.0)
print("Total Number of Ca Between -5.0 - 5.0   :", plane, " across ", data['timesteps'], "timesteps")
print("Average Number of Ca Between -5.0 - 5.0 :", plane / data['timesteps'])

plane = f_density.one_dimensional_density_sb(ul=5.0, ll=-5.0)
print("Total Number of F Between -5.0 - 5.0   :", plane, " across ", data['timesteps'], "timesteps")
print("Average Number of F Between -5.0 - 5.0 :", plane / data['timesteps'])