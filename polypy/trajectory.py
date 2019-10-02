import numpy as np 

class PolyTrajectory():
    """The PolyTrajectory class is a collection of functions for manipulation
    of the data generated from the HISTORY file.

    Parameters
    ----------
    traj : dictionary
        dictionary containing data from the entire history file.
    """
    def __init__(self, traj):
        self.traj = traj


    def get_file_type(self):
        """Return the filetype
        
        Returns
        -------
        str
            File type
        """
        return self.traj['trajectory_type']

    
    def lattice_vectors(self):
        """Returns the unique elements of the lattice vectors.
        Only valid for an orthogonal system.
        
        Returns
        -------
        iconfig_lvs : array like
            array of a,b,c lattice vectors.
        """
        iconfig_lvs = np.array([])

        for i in range(self.traj['numconfigs']):
            iconfig_lvs = np.append(iconfig_lvs, self.traj[i]['lvs'])
        iconfig_lvs = np.split(iconfig_lvs, self.get_nconfigs())
        for i in range(len(iconfig_lvs)):
            iconfig_lvs[i] = np.reshape(iconfig_lvs[i], (3, 3))
        iconfig_lvs = np.asarray(iconfig_lvs)
        return iconfig_lvs
                

    def get_nconfigs(self):
        """Returns the total number of configs found in the HISTORY file.

        Returns
        -------
        int
            Total number of configs in the history file.
        """
        return int(self.traj['numconfigs']) 
   

    def get_nconfigs(self):
        """Returns the total number of configs found in the HISTORY file.

        Returns
        -------
        int
            Total number of configs in the history file.
        """
        return int(self.traj['numconfigs']) 
   

  #  def timestep(self):
  #      """Calculates the simulation timestep.#

#        Returns
 #       -------
  #      float
   #         Timestep - Time between records.
    #    """
     #   timestep = float(self.traj[0]['timestep'][4])
      #  step_length = float(self.traj[1]['timestep'][0]) - float(self.traj[0]['timestep'][0])
       # return (timestep * step_length)

 #       if self.traj['trajectory_type'] != "DLPOLY":
  #          print("Monte Carlo does not have a timestep")
   #     else:       
    #        timestep = float(self.traj[0]['timestep'][4])
     #       step_length = float(self.traj[1]['timestep'][0]) - float(self.traj[0]['timestep'][0])
      #      return (timestep * step_length)

    def get_title(self):
        """Returns the title of the HISTORY file.

        Returns
        -------
        list
            List of elements in the HISTORY file title.
        """
        if self.traj['trajectory_type'] == "DLPOLY":
            return self.traj['title']
        else:
            print("No title in a DLMONTE file")


    def get_style(self):
        """Returns the HISTORY file style.

        Returns
        -------
        list
            List of elements on the style line of the history file.
        """
        if self.traj['trajectory_type'] == "DLPOLY":
            return self.traj['style']
        else:
            print("No style in a DLMONTE file")


    def get_total_natoms(self):
        """Returns the total number of atoms in the simulation

        Returns
        -------
        int
            Total number of atoms in the simulation
        """
        if self.traj['trajectory_type'] == "DLPOLY":
            return self.get_dlpoly_total_natoms()
        elif self.traj['trajectory_type'] == "DLMONTE":
            return self.get_dlmonte_total_natoms()


    def atom_coordinates(self, atom_label):
        """Isolate the coordinates for a specific atom.

        Parameters
        ----------
        atom_label : str
            Label of the desired atom.

        Returns
        -------
        atom_coords : array like
            array of x,y,z coordinates for every desired atom in the simulation.
        """
        if self.traj['trajectory_type'] == "DLPOLY":
            return self.get_dlpoly_atom_coordinates(atom_label)
        elif self.traj['trajectory_type'] == "DLMONTE":
            return self.get_dlmonte_atom_coordinates(atom_label)


    def get_natoms(self, atom_label):
        """Returns the total number of a specific atom in the simulation.

        Parameters
        ----------
        atom_label : str
           Label of the desired atom.       
        
        Returns
        -------
        int
            Total number of desired atoms in the simulation
        """
        if self.traj['trajectory_type'] == "DLPOLY":
            return self.get_dlpoly_natom(atom_label)
        elif self.traj['trajectory_type'] == "DLMONTE":
            return self.get_dlmonte_natom(atom_label)


    def get_dlmonte_total_natoms(self):
        natoms = np.array([])
        for i in range(self.traj['numconfigs']):
            config_natoms = np.array([])
            for j in range(self.traj[i]['nummols']):
                config_natoms = np.append(config_natoms, self.traj[i]['mols'][j]['numatoms'])
                natoms = np.append(natoms, np.sum(config_natoms))
        if np.mean(natoms) == natoms[0]:
           print("Atom number does not change")
           return np.mean(natoms)
        else:
            print("Atom number changes")
            return natoms


    def get_dlpoly_total_natoms(self):
        return self.traj[0]['atoms']['numatoms']


    def get_dlpoly_atom_coordinates(self, atom_label):
        atom_coords = []
        for i in range(self.traj['numconfigs']):
            for j in range(self.traj[i]['atoms']['numatoms']):
                if self.traj[i]['atoms'][j]['label'] == atom_label:
                    atom_coords.append(self.traj[i]['atoms'][j]['coor'])
        atom_coords = np.asarray(atom_coords)
        atom_coords = np.split(atom_coords, self.traj['numconfigs'])
        return atom_coords


    def get_dlmonte_atom_coordinates(self, atom_label):
        atom_coords = []
        for i in range(self.traj['numconfigs']):
            for j in range(self.traj[i]['nummols']):
                for k in range(self.traj[i]['mols'][j]['maxatoms']):
                    if self.traj[i]['mols'][j]['atoms'][k]['label'] == atom_label:
                        atom_coords.append(self.traj[i]['mols'][j]['atoms'][k]['coor'])
        atom_coords = np.asarray(atom_coords)
        return atom_coords


    def get_dlpoly_natom(self, atom_label):
        atom_count = 0
        for i in range(self.traj['numconfigs']):
            for j in range(self.traj[i]['atoms']['numatoms']):
                if self.traj[i]['atoms'][j]['label'] == atom_label:
                    atom_count = atom_count + 1
        return int(atom_count / self.traj['numconfigs'])


    def get_dlmonte_natom(self, atom_label):
        atom_count = 0
        for i in range(self.traj['numconfigs']):
            for j in range(self.traj[i]['nummols']):
                for k in range(self.traj[i]['mols'][j]['maxatoms']):
                    if self.traj[i]['mols'][j]['atoms'][k]['label'] == atom_label:
                       atom_count = atom_count + 1
        return int(atom_count / self.traj['numconfigs'])
