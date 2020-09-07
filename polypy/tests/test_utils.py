import numpy as np
from polypy import utils as ut
from polypy import trajectory 
import unittest
from numpy.testing import assert_almost_equal
import os


test_history = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/HISTORY1')


class TestUtils(unittest.TestCase):

    def test_system_volume(self):
        traj = dlppy.read_trajectory(test_history)
        traj_ob = trajectory.PolyTrajectory(traj)
        volume = ut.volume

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
                                 np.array([5, 10, 15, 20, 25,
                                          30, 35, 40, 45]))[0]
        expected = 5.0
        assert a == expected

    def test_pbc(self):
        a, b = ut.pbc(10, 9, 20)
        c, d = ut.pbc(1, 40, 20)
        expected_a = False
        expected_b = 10
        expected_c = True
        expected_d = 41
        assert a == expected_a
        assert b == expected_b
        assert c == expected_c
        assert d == expected_d

    def test_bin_choose(self):
        a = ut.bin_choose(4, 2)
        expected = 1
        assert a == expected

    def test_one_dimensional_charge_density(self):

        atoms_coords = [np.array([1, 3, 5, 7, 9]),
                        np.array([2, 4, 6, 8, 10])]
        atom_charges = [1.0, -1.0]
        charge_density = ut.one_dimensional_charge_density(atoms_coords,
                                                           atom_charges,
                                                           1.0)
        assert_almost_equal(charge_density, np.array([-1.0,
                                                      -1.0,
                                                      -1.0,
                                                      -1.0,
                                                      -1.0]))

    def test_two_dimensional_charge_density(self):

        atoms_coords = [np.array([[1, 3, 5],
                                  [7, 9, 11]]),
                        np.array([[2, 4, 6],
                                  [8, 10, 12]])]
        atom_charges = [1.0, -1.0]
        charge_density = ut.two_dimensional_charge_density(atoms_coords,
                                                           atom_charges,
                                                           1.0)
        assert_almost_equal(charge_density, np.array([[-1.0, -1.0, -1.0],
                                                      [-1.0, -1.0, -1.0]]))

    def test_poisson_solver(self):
        dx, e_field, potential = ut.poisson_solver(np.array([1, 2, 3, 4, 5]),
                                                   np.array([-1.0,
                                                             -1.0,
                                                             -1.0,
                                                             -1.0,
                                                             -1.0]),
                                                   1)
        predicted_e_field = np.array([2.87995168e+01,
                                      1.43997584e+01,
                                      3.55271368e-15,
                                      -1.43997584e+01,
                                      -2.87995168e+01])
        predicted_potential = np.array([-0.00000000e+00,
                                        -2.15996376e+01,
                                        -2.87995168e+01,
                                        -2.15996376e+01,
                                        -1.06581410e-14])
        assert_almost_equal(e_field, predicted_e_field)
        assert_almost_equal(potential, predicted_potential)

    def test_smooth_msd_data(self):
        x_data = np.array([1, 2, 3, 4, 5, 1, 2, 3, 4, 5])
        y_data = np.array([5, 4, 3, 2, 1, 5, 4, 3, 2, 1])
        x, y = ut.smooth_msd_data(x_data, y_data)
        predicted_x = np.array([1, 2, 3, 4, 5])
        predicted_y = np.array([5, 4, 3, 2, 1])
        assert_almost_equal(x, predicted_x)
        assert_almost_equal(y, predicted_y)
