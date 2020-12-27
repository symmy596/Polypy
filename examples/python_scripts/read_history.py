from polypy import read as rd

history = rd.History("../example_data/HISTORY_CaF2", ["CA", "F"])

print(history.trajectory.fractional_trajectory)

print(history.trajectory.timesteps)

print(history.trajectory.atoms_in_history)

print(history.trajectory.total_atoms)