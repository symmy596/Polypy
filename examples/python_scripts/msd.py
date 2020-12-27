from polypy import read as rd
from polypy.msd import MSD 
from polypy.msd import RegionalMSD 
from polypy import analysis
from polypy import utils as ut
from polypy import plotting
import numpy as np
import matplotlib.pyplot as plt

history_caf2 = rd.History("../example_data/HISTORY_CaF2", ["F"])

f_msd = MSD(history_caf2.trajectory, sweeps=2)

output = f_msd.msd()

ax = plotting.msd_plot(output)

plt.savefig("msd.png", dpi=600)

f_msd = MSD(history_caf2.trajectory, sweeps=10)

output = f_msd.msd()

ax = plotting.msd_plot(output)
plt.savefig("msd_2.png", dpi=600)

print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())

volume, step = analysis.system_volume(history_caf2.trajectory)
average_volume = np.mean(volume[:50])

sigma = analysis.conductivity(history_caf2.trajectory.total_atoms, 
                        average_volume, 
                        output.xyz_diffusion_coefficient(), 
                        1500, 1)

print("Ionic Conductivity :", sigma)

print("Resistivity :", (1 / sigma)) 
