from polypy import read as rd
from polypy import density as Dens
from polypy import utils as ut
from polypy import write as wr
import numpy as np


data = rd.read_history("../example_data/HISTORY", ["CA", "F"])

ca_density = Dens.Density(data, atom_type="CA")
f_density = Dens.Density(data, atom_type="F")

cx, cy = ca_density.one_dimensional_density(Bin=0.1, direction="x")
fx, fy = f_density.one_dimensional_density(Bin=0.1, direction="x")

bin_volume = 0.1 * np.mean(data['lv'][:,1] * np.mean(data['lv'][:,2]))

charge_density = ut.one_dimensional_charge_density([fy, cy], [-1.0, 2.0], bin_volume)

wr.one_dimensional_charge_density_plot(fx, charge_density)