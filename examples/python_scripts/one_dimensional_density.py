from polypy.read import History
from polypy.read import Archive
from polypy.density import Density
from polypy import analysis
from polypy import utils as ut
from polypy import plotting

import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

history = History("../example_data/HISTORY_GB", ["CE", "O"])

ce_density = Density(history.trajectory, atom="CE", histogram_size=0.1)
o_density = Density(history.trajectory, atom="O", histogram_size=0.1)

cx, cy, c_volume = ce_density.one_dimensional_density(direction="z")
ox, oy, o_volume = o_density.one_dimensional_density(direction="z")

ax = plotting.one_dimensional_density_plot([cx, ox], [cy, oy], ["Ce", "O"])
ax.set_xlim(42, 82)
plt.savefig("one_dimensional_density.png", dpi=600)

charge = analysis.OneDimensionalChargeDensity(ox, [oy, cy], [-2.0, 4.0], c_volume, history.trajectory.timesteps)

dx, charge_density = charge.calculate_charge_density()

ax = plotting.one_dimensional_charge_density_plot(dx, charge_density)
ax.set_xlim(42, 82)

plt.savefig("one_dimensional_charge_density.png", dpi=600)


dx, electric_field = charge.calculate_electric_field()

ax = plotting.electric_field_plot(dx, electric_field)
ax.set_xlim(42, 82)
plt.savefig("electric_field.png", dpi=600)


dx, electrostatic_potential = charge.calculate_electrostatic_potential()

ax = plotting.electrostatic_potential_plot(dx, electrostatic_potential)
ax.set_xlim(42, 82)

plt.savefig("electrostatic_potential.png", dpi=600)
