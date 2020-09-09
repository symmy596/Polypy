import numpy as np
import os
from polypy import read as rd
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')
test_config = os.path.join(os.path.dirname(__file__), 'CONFIG')
test_archive = os.path.join(os.path.dirname(__file__), 'ARCHIVE')


class TestHistory(unittest.TestCase):

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
        data = rd.History(test_history, ["CA"])
        assert_almost_equal(expected_traj, data.trajectory.cartesian_trajectory)
        assert data.trajectory.total_atoms == 1
        assert data.trajectory.timesteps == 10

class TestConfig(unittest.TestCase):

    def test_read_config(self):
        expected_traj = np.array([[1.0, 1.0, 1.0]])
        data = rd.Config(test_config, ["CA"])
        assert_almost_equal(expected_traj, data.trajectory.cartesian_trajectory)
        assert data.trajectory.total_atoms == 1
        assert data.trajectory.timesteps == 1


class TestArchive(unittest.TestCase):

    def test_read_archive(self):
        expected_traj = np.array([[1.0, 1.0, 1.0],
                                  [2.0, 2.0, 2.0],
                                  [3.0, 3.0, 3.0]])
        data = rd.Archive(test_archive, ["Ca"])
        assert_almost_equal(expected_traj, data.trajectory.cartesian_trajectory)

class TestTrajectory(unittest.TestCase):

    def test_get_config(self):
        data = rd.History(test_history, ["CA", "F"])
        expected_config = data.trajectory.get_config(0)
        predicted_config = np.array([[1.00, 1.00, 1.00],
                                     [10.00, 10.00, 10.00]])
        assert_almost_equal(expected_config.cartesian_trajectory, predicted_config)

    def test_get_atom(self):
        data = rd.History(test_history, ["CA", "F"])
        expected_ca = data.trajectory.get_atom("CA")
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
        assert_almost_equal(predicted_ca, expected_ca.cartesian_trajectory)
