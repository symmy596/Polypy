import numpy as np
import numpy.testing as npt
import os
from polypy import read_dl_poly as dlppy
import unittest
from numpy.testing import assert_almost_equal
 
test_history = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/HISTORY1')
test_config1 = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/CONFIG1')
test_config2 = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/CONFIG2')
test_config3 = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/CONFIG3')

class TestReadDlpoly(unittest.TestCase):

    def test_open_config_numatoms(self):
        config = dlppy.open_config(test_config1)
        expect_numatoms = 3
        assert config['atoms']['numatoms'] == expect_numatoms

    def test_open_config_atomlabels(self):
        config = dlppy.open_config(test_config1)
        expect_atom_labels = ["A", "B", "C"]
        for iatom in range ( config['atoms']['numatoms'] ):
            assert config['atoms'][iatom]['label'] == expect_atom_labels[iatom]

    def test_open_config_atomcoors(self):
        config = dlppy.open_config(test_config1)
        expect_atom_coors = [ [0.1,0.2,0.3], 
                            [0.5, 0.6, 0.7], 
                            [1.0, 1.5, 2.0] ]
        for iatom in range ( config['atoms']['numatoms'] ):
            npt.assert_array_almost_equal( 
                config['atoms'][iatom]['coor'], 
                expect_atom_coors[iatom] )

    def test_open_config_atomcoors2(self):
        config = dlppy.open_config(test_config2)
        expect_atom_coors = [ [0.1,0.2,0.3], 
                            [0.5, 0.6, 0.7], 
                            [1.0, 1.5, 2.0] ]
        for iatom in range ( config['atoms']['numatoms'] ):
            npt.assert_array_almost_equal( 
                config['atoms'][iatom]['coor'], 
                expect_atom_coors[iatom] )

    def test_open_config_atomvels2(self):
        config = dlppy.open_config(test_config2)
        expect_atom_vels = [ [0.01,0.02,0.03], 
                            [0.05, 0.06, 0.07], 
                            [1.01, 1.51, 2.01] ]
        for iatom in range ( config['atoms']['numatoms'] ):
            npt.assert_array_almost_equal( 
                config['atoms'][iatom]['vel'], 
                expect_atom_vels[iatom] )

    def test_open_config_atomcoors3(self):
        config = dlppy.open_config(test_config3)
        expect_atom_coors = [ [0.1,0.2,0.3], 
                            [0.5, 0.6, 0.7], 
                            [1.0, 1.5, 2.0] ]
        for iatom in range ( config['atoms']['numatoms'] ):
            npt.assert_array_almost_equal( 
                config['atoms'][iatom]['coor'], 
                expect_atom_coors[iatom] )

    def test_open_config_atomvels3(self):
        config = dlppy.open_config(test_config3)
        expect_atom_vels = [ [0.01,0.02,0.03], 
                            [0.05, 0.06, 0.07], 
                            [1.01, 1.51, 2.01] ]
        for iatom in range ( config['atoms']['numatoms'] ):
            npt.assert_array_almost_equal( 
                config['atoms'][iatom]['vel'], 
                expect_atom_vels[iatom] )

    def test_open_config_atomforces3(self):
        config = dlppy.open_config(test_config3)
        expect_atom_forces = [ [0.001,0.002,0.003], 
                            [0.005, 0.006, 0.007], 
                            [1.001, 1.501, 2.001] ]
        for iatom in range ( config['atoms']['numatoms'] ):
            npt.assert_array_almost_equal( 
                config['atoms'][iatom]['force'], 
                expect_atom_forces[iatom] )


    def test_read_archive_numconfigs(self):
        traj = dlppy.read_trajectory(test_history)
        expect_configs = 3
        assert traj['numconfigs'] == expect_configs

    def test_read_archive_atomlabels(self):
        traj = dlppy.read_trajectory(test_history)
        expect_labels = [ [ "A", "B", "C" ] , 
                        [ "D", "E", "F" ] ,
                        [ "G", "H", "I" ] ]
        for iconfig in range( traj['numconfigs']):
            this_config = traj[iconfig]
            for iatom in range(this_config['atoms']['numatoms']):
                    this_atom = this_config['atoms'][iatom]
                    assert ( this_atom['label'] == 
                        expect_labels[iconfig][iatom] )
