from polypy import analysis
from polypy import plotting
from polypy import read as rd
import matplotlib.pyplot as plt

history = rd.Archive("../example_data/Archive_LLZO", ["ZR"])

volume, step = analysis.system_volume(history.trajectory)
ax = plotting.volume_plot(step, volume)
plt.savefig("archive_volume.png", dpi=600)
