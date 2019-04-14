import numpy as np
import numpy.testing as npt
import os
from polypy import read_dl_monte as dlmpy
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/ARCHIVE1')
test_config1 = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/CONFIG1')

class TestReadDlpoly(unittest.TestCase):

    def test_open_config_matchnummols(self):
        '''Read entire configuration from CONFIG
        Check that molecules numbers match'''

        config_file = test_data_dir+'sample_configs/CONFIG1'

        config = dlmpy.open_config(config_file)

        assert config['nummols'] == len(config['mols'])

    def test_open_config_moleculenames(self):
        '''Read entire configuration from CONFIG
        Check that molecules names match'''

        config_file = test_data_dir+'sample_configs/CONFIG1'

        config = dlmpy.open_config(config_file)

        expect_mol_names = [ 'bob', 'charlie' ]

        for imol in range( config['nummols']):
            assert config['mols'][imol]['name']== expect_mol_names[imol]

    def test_open_config_atomnums(self):
        '''Read entire configuration from CONFIG
        Check that atom numberss match'''

        config_file = test_data_dir+'sample_configs/CONFIG1'

        config = dlmpy.open_config(config_file)

        expect_numatoms = [1, 2]

        for imol in range( config['nummols']):
            assert config['mols'][imol]['numatoms'] == expect_numatoms[imol]

    def test_open_config_atomlabels(self):
        '''Read entire configuration from CONFIG
        Check that atom labels match'''

        config_file = test_data_dir+'sample_configs/CONFIG1'

        config = dlmpy.open_config(config_file)

        expect_atom_labels = {0: ["A"], 1: ["B", "C"] }

        for imol in range( config['nummols']):
            for iatom in range ( config['mols'][imol]['numatoms'] ):
                assert config['mols'][imol]['atoms'][iatom]['label'] == expect_atom_labels[imol][iatom]

    def test_open_config_atomcoors(self):
        '''Read entire configuration from CONFIG
        Check that atom coors match'''

        config_file = test_data_dir+'sample_configs/CONFIG1'

        config = dlmpy.open_config(config_file)

        expect_atom_coors = {0: [ np.array([0.1,0.2,0.3]) ], 
                            1: np.array([ [0.5, 0.6, 0.7], [1.0, 1.5, 2.0] ]) }

        for imol in range( config['nummols']):
            for iatom in range ( config['mols'][imol]['numatoms'] ):
                npt.assert_array_almost_equal( 
                config['mols'][imol]['atoms'][iatom]['coor'], 
                expect_atom_coors[imol][iatom] )

    def test_read_archive_numconfigs(self):
        '''Read entire trajectory from ARCHIVE
        Check that number of configurations match'''

        archive_file = test_data_dir+'sample_configs/ARCHIVE1'

        traj = dlmpy.read_trajectory(archive_file)

        expect_configs = 3

        assert traj['numconfigs'] == expect_configs

    def test_read_archive_atomlabels(self):
        '''Read entire trajectory from ARCHIVE
        Check that labels of atoms in configurations match'''

        archive_file = test_data_dir+'sample_configs/ARCHIVE1'

        traj = dlmpy.read_trajectory(archive_file)

        expect_labels = [ [ ["A"], ["B", "C"] ] , 
                        [ ["D"], ["E", "F"] ] ,
                        [ ["G"], ["H", "I"] ] ]

        for iconfig in range( traj['numconfigs']):
            this_config = traj[iconfig]
            for imol in range ( this_config['nummols'] ):
                this_mol = this_config['mols'][imol]
                for iatom in range(this_mol['numatoms']):
                    this_atom = this_mol['atoms'][iatom]
                    assert ( this_atom['label'] == 
                        expect_labels[iconfig][imol][iatom] )
