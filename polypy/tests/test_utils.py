import numpy as np
from polypy import Utils as ut
import unittest
from numpy.testing import assert_almost_equal
import os


test_volume = {'label': np.array(["H", "H", "H", "H", "H", "H", "H", "H", "H", "H",]), 
        'trajectories': np.array([[15.00,	15.00,	15.00],
                                  [16.00,	16.00,	16.00],
                                  [17.00,	17.00,	17.00],
                                  [18.00,	18.00,	18.00],
                                  [19.00,	19.00,	19.00],
                                  [20.00,	20.00,	20.00],
                                  [21.00,	21.00,	21.00],
                                  [22.00,	22.00,	22.00],
                                  [23.00,	23.00,	23.00],
                                  [24.00,	24.00,	24.00]]),
        'lv': np.array([[15.00,	15.00,	15.00],
                        [16.00,	16.00,	16.00],
                        [17.00,	17.00,	17.00],
                        [18.00,	18.00,	18.00],
                        [19.00,	19.00,	19.00],
                        [20.00,	20.00,	20.00],
                        [21.00,	21.00,	21.00],
                        [22.00,	22.00,	22.00],
                        [23.00,	23.00,	23.00],
                        [24.00,	24.00,	24.00]]),
        'timesteps': 10,
         'natoms': 10}

class TestUtils(unittest.TestCase):

    def test_system_volume(self):

        expected_vol = np.array([3375.000, 4096.000, 4913.000, 5832.000, 6859.000, 8000.000, 9261.000, 10648.000, 12167.000, 13824.000])
        expected_time = np.arange(10) * 10

        a, b = ut.system_volume(test_volume, 10)
        assert_almost_equal(expected_vol, a)
        assert_almost_equal(expected_time, b)
    
    def test_conductivity(self):
        a = ut.conductivity(1, 2, 3, 4)
        a = round(a, 0)
        expected = 697219
        assert_almost_equal(a, expected)
    
    def test_three_d_diffusion_coefficient(self):
        a = ut.three_d_diffusion_coefficient(10)
        expected = 16.666666666
        assert_almost_equal(a, expected)
    
    def test_one_d_diffusion_coefficient(self):
        a = ut.one_d_diffusion_coefficient(10)
        expected = 50.0
        assert a == expected

    def test_linear_regression(self):
        a = ut.linear_regression(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                 np.array([5, 10, 15, 20, 25, 30, 35, 40, 45]))
        expected = 5.0
        assert a == expected

    def test_charge_density(self):
        pass
    
    def test_poisson_solver(self):
        pass
