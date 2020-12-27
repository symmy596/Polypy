
from polypy import read as rd

archive = rd.Archive("../example_data/ARCHIVE_Short", ["AL"])

print(archive.trajectory.fractional_trajectory)

print(archive.trajectory.timesteps)

print(archive.trajectory.atoms_in_history)

print(archive.trajectory.total_atoms)
