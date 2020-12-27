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

cx_2d, cy_2d, cz_2d, c_volume = ce_density.two_dimensional_density(direction="x")
ox_2d, oy_2d, oz_2d, o_volume = o_density.two_dimensional_density(direction="x")

fig, ax = plotting.two_dimensional_density_plot(cx_2d, cy_2d, cz_2d, colorbar=False, palette="Greys", log=True)
ax.set_xlim(42, 82)
ax.axis('off')
plt.savefig("ce_two_dimensional_density.png", dpi=600)

fig, ax = plotting.two_dimensional_density_plot(ox_2d, oy_2d, oz_2d, colorbar=False, palette="Greys", log=True)
ax.set_xlim(42, 82)
ax.axis('off')
plt.savefig("o_two_dimensional_density.png", dpi=600)


charge_density = analysis.two_dimensional_charge_density([oz_2d, cz_2d], [-2.0, 4.0], o_volume, history.trajectory.timesteps)
fig, ax = plotting.two_dimensional_charge_density_plot(ox_2d, oy_2d, charge_density, palette='bwr')
ax.set_xlim(42, 82)
plt.savefig("two_dimensional_charge_density.png", dpi=600)



fig, ax = plotting.combined_density_plot(cx_2d, cy_2d, cz_2d, palette="Oranges", linecolor="orange", log=True)
for axes in ax:
    axes.set_xlim(42, 82)
plt.savefig("ce_combined_two_dimensional_density.png", dpi=600)

fig, ax = plotting.combined_density_plot(ox_2d, oy_2d, oz_2d, palette="Blues", linecolor="blue", log=True)
for axes in ax:
    axes.set_xlim(42, 82)
plt.savefig("o_combined_two_dimensional_density.png", dpi=600)


fig, ax = plotting.two_dimensional_density_plot_multiple_species([cx_2d, ox_2d], [cy_2d, oy_2d], 
                                                                 [cz_2d, oz_2d], ["Blues", "Oranges"], 
                                                                 log=True)
ax.set_xlim(42, 82)
plt.savefig("combined_two_dimensional_density.png", dpi=600)


fig, ax = plotting.combined_density_plot_multiple_species(x_list=[cx_2d, ox_2d], 
                                                          y_list=[cy_2d, oy_2d], 
                                                          z_list=[cz_2d, oz_2d], 
                                                          palette_list=["Blues", "Oranges"], 
                                                          label_list=['Ce', 'O'], 
                                                          color_list=["blue", "orange"],
                                                          log=True)
for axes in ax:
    axes.set_xlim(42, 82)
plt.savefig("combined_two_dimensional_density_mp.png", dpi=600)
