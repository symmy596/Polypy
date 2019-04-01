import numpy as np
import os
from polypy import read as rd
from polypy import density as dens
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')
test_config = os.path.join(os.path.dirname(__file__), 'CONFIG')
test_archive = os.path.join(os.path.dirname(__file__), 'ARCHIVE')
expected_z = os.path.join(os.path.dirname(__file__), 'Expected_Z')


class TestDensity(unittest.TestCase):

    def test_one_dimensional_density(self):
        data = rd.read_history(test_history, ["CA"])
        test_density = dens.Density(data)
        xx, yx = test_density.one_dimensional_density(histogram_width=1.00,
                                                      direction="x")
        xy, yy = test_density.one_dimensional_density(histogram_width=1.00,
                                                      direction="y")
        xz, yz = test_density.one_dimensional_density(histogram_width=1.00,
                                                      direction="z")

        predicted_x = np.array([-5.0, -4.0, -3.0, -2.0, -1.0,
                                0.0, 1.0, 2.0, 3.0, 4.0])
        predicted_y = np.zeros(10) + 1
        assert_almost_equal(xx, predicted_x)
        assert_almost_equal(yx, predicted_y)
        assert_almost_equal(xy, predicted_x)
        assert_almost_equal(yy, predicted_y)
        assert_almost_equal(xz, predicted_x)
        assert_almost_equal(yz, predicted_y)

    def test_one_dimensional_density_sb(self):
        data = rd.read_history(test_history, ["CA"])
        test_density = dens.Density(data)
        plane = test_density.one_dimensional_density_sb(ll=-2.0, ul=2.0)
        predicted_plane = 3.0
        assert plane == predicted_plane

    def test_two_dimensional_density(self):
        data = rd.read_history(test_history, ["CA"])
        test_density = dens.Density(data)
        x, y, z = test_density.two_dimensional_density(box=1.0)
        predicted_x = np.array([-5.0, -4.0, -3.0, -2.0,
                                -1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
        predicted_y = np.array([0.0, 1.0, 2.0, 3.0, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0])
        predicted_z = np.genfromtxt(expected_z, delimiter=",", dtype="float")
        assert_almost_equal(x, predicted_x)
        assert_almost_equal(y, predicted_y)
        assert_almost_equal(z, predicted_z)

    def test_one_and_two_dimension_overlay(self):
        data = rd.read_history(test_history, ["CA"])
        test_density = dens.Density(data)
        x, y, z, y2 = test_density.one_and_two_dimension_overlay(box=1.0)
        predicted_x = np.array([-5.0, -4.0, -3.0, -2.0,
                                -1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
        predicted_y = np.array([0.0, 1.0, 2.0, 3.0, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0])
        predicted_z = np.genfromtxt(expected_z, delimiter=",", dtype="float")
        predicted_y2 = np.zeros(10) + 1.001
        assert_almost_equal(x, predicted_x)
        assert_almost_equal(y, predicted_y)
        assert_almost_equal(z, predicted_z)
        assert_almost_equal(y2, predicted_y2)
