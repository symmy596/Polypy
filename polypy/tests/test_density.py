import numpy as np
import os
from polypy import read
from polypy.density import Density
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')
expected_z = os.path.join(os.path.dirname(__file__), 'Expected_Z')


class TestDensity(unittest.TestCase):

    def test_find_limits(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        assert test_density.x_lim == 10
        assert test_density.y_lim == 10
        assert test_density.z_lim == 10
        assert test_density.x.size == 10

    def test_build_map(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        assert np.sum(test_density.coords_map) == 10.0
        test_density.build_map()
        assert np.sum(test_density.coords_map) == 20.0

    def test_update_map(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        test_density.update_map([0.1, 0.1, 0.1])
        assert test_density.coords_map[1, 1, 1] == 2.0

    def test_one_dimensional_density(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        xx, yx, vol = test_density.one_dimensional_density(direction="x")

        predicted_x = np.array([0.0, 1.0, 2.0, 3.0, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0])
        predicted_y = np.zeros(10) + 1
        assert_almost_equal(xx, predicted_x)
        assert_almost_equal(yx, predicted_y)

    def test_two_dimensional_density(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        x, y, z, vol = test_density.two_dimensional_density()
        predicted_x = np.array([0.0, 1.0, 2.0, 3.0, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0])
        predicted_y = np.array([0.0, 1.0, 2.0, 3.0, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0])
        predicted_z = np.genfromtxt(expected_z, delimiter=",", dtype="float")
        assert_almost_equal(x, predicted_x)
        assert_almost_equal(y, predicted_y)
        assert_almost_equal(z, predicted_z)
