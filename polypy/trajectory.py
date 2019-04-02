import numpy as np 

class Trajectory():

    def __init__(self, traj):
        self.traj = traj


    def lattice_vectors(self):
        iconfig_lvs = np.array([])
        for i in range(self.traj['numconfigs']):
            iconfig_lvs = np.append(iconfig_lvs, np.sum(self.traj[i]['lvs'], axis=1))
        iconfig_lvs = np.reshape(iconfig_lvs, (self.traj['numconfigs'], 3))
        return iconfig_lvs


    def timestep(self):
        timestep = float(self.traj[0]['timestep'][4])
        step_length = float(self.traj[1]['timestep'][0]) - float(self.traj[0]['timestep'][0])
        return (timestep * step_length)


    def get_nconfigs(self):
        return self.traj['numconfigs'] 

    def get_title(self):
        return self.traj['title']


    def get_style(self):
        return self.traj['style']


    def get_total_natoms(self):
        natoms = np.array([])
        for i in range(self.traj['numconfigs']):
            natoms = np.append(natoms, self.traj[i]['atoms']['numatoms'])
        if natoms[0] == np.mean(natoms):
            return np.mean(natoms)
        else:
            print("The number of atoms changes between configs")
            return natoms

    
    def atom_coordinates(self, atom_label):
        atom_coords = []
        for i in range(self.traj['numconfigs']):
            for j in range(self.traj[i]['atoms']['numatoms']):
                if self.traj[i]['atoms'][j]['label'] == atom_label:
                    atom_coords.append(self.traj[i]['atoms'][j]['coor'])
        atom_coords = np.asarray(atom_coords)
        return atom_coords


    def get_natoms(self, atom_label):
        atom_count = 0
        for i in range(self.traj['numconfigs']):
            for j in range(self.traj[i]['atoms']['numatoms']):
                if self.traj[i]['atoms'][j]['label'] == atom_label:
                    atom_count = atom_count + 1
        return (atom_count / self.traj['numconfigs'])
    
    
