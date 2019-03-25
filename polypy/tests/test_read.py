import numpy as np
import os
from polypy import read as rd
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')
test_config = os.path.join(os.path.dirname(__file__), 'CONFIG')
test_archive = os.path.join(os.path.dirname(__file__), 'ARCHIVE')

class TestUtils(unittest.TestCase):

    def test_read_history(self):
        expected_traj = np.array([[1.0, 1.0, 1.0],
                                  [2.0, 2.0, 2.0],
                                  [3.0, 3.0, 3.0],
                                  [4.0, 4.0, 4.0],
                                  [0.0, 0.0, 0.0],
                                  [-1.0, -1.0, -1.0],
                                  [-2.0, -2.0, -2.0],
                                  [-3.0, -3.0, -3.0],
                                  [-4.0, -4.0, -4.0],
                                  [-5.0, -5.0, -5.0]])
        data = rd.read_history(test_history, ["CA"])
        assert_almost_equal(expected_traj, data['trajectories'])
        assert data['natoms'] == 1
        assert data['timesteps'] == 10


    def test_read_config(self):
        expected_traj = np.array([[1.0, 1.0, 1.0]])
        data = rd.read_config(test_config, ["CA"])
        assert_almost_equal(expected_traj, data['trajectories'])
        assert data['natoms'] == 1
        assert data['timesteps'] == 1

    def test_read_archive(self):
        expected_traj = np.array([[1.0, 1.0, 1.0, 0],
                                  [2.0, 2.0, 2.0, 0],
                                  [3.0, 3.0, 3.0, 0]])
        data = rd.read_archive(test_archive, ["Ca"])
        assert_almost_equal(expected_traj, data['trajectories'])
        assert data['natoms'] == 1
        assert data['timesteps'] == 3

    def test_get_config(self):
        data = rd.read_history(test_history, ["CA", "F"])
        expected_config = rd.get_config(data, 0)
        predicted_config = np.array([[1.00, 1.00, 1.00],
                                      [10.00, 10.00, 10.00]])
        assert_almost_equal(expected_config, predicted_config)


    def test_get_atom(self):        
        data = rd.read_history(test_history, ["CA", "F"])
        expected_ca = rd.get_atom(data, "CA")
        predicted_ca = np.array([[1.0, 1.0, 1.0],
                                  [2.0, 2.0, 2.0],
                                  [3.0, 3.0, 3.0],
                                  [4.0, 4.0, 4.0],
                                  [0.0, 0.0, 0.0],
                                  [-1.0, -1.0, -1.0],
                                  [-2.0, -2.0, -2.0],
                                  [-3.0, -3.0, -3.0],
                                  [-4.0, -4.0, -4.0],
                                  [-5.0, -5.0, -5.0]])
        assert_almost_equal(predicted_ca, expected_ca['trajectories'])        
