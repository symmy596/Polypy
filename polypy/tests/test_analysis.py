import numpy as np
import os
import unittest
from polypy import analysis
from polypy.density import Density
from polypy import read
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')


class TestAnalysis(unittest.TestCase):

    def test_conductivity(self):
        a = analysis.conductivity(1, 2, 3, 4, 1)
        a = round(a, 0)
        expected = 7
        assert_almost_equal(a, expected)

    def test_system_volume(self):
        data = read.History(test_history, ['CA'])
        vol = analysis.system_volume(data.trajectory)[0]
        expected = np.zeros(10) + 1000
        assert_almost_equal(vol, expected)

    def test_two_dimensional_charge_density(self):
        atoms_coords = [np.array([[1, 3, 5],
                                  [7, 9, 11]]),
                        np.array([[2, 4, 6],
                                  [8, 10, 12]])]
        atom_charges = [1.0, -1.0]
        charge_density = analysis.two_dimensional_charge_density(atoms_coords,
                                                                 atom_charges,
                                                                 1.0, 1.0)
        assert_almost_equal(charge_density, np.array([[-1.0, -1.0, -1.0],
                                                      [-1.0, -1.0, -1.0]]))


class TestOneDimensionalChargeDensity(unittest.TestCase):

    def test_calculate_charge_density(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        xx, yx, vol = test_density.one_dimensional_density(direction="x")
        charge = analysis.OneDimensionalChargeDensity(xx, [yx], [2.0], vol, 1)
        dx, charge_density = charge.calculate_charge_density()
        expected = np.zeros(10) + 0.2
        assert_almost_equal(charge_density, expected)

    def test_calculate_electric_field(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        xx, yx, vol = test_density.one_dimensional_density(direction="x")
        charge = analysis.OneDimensionalChargeDensity(xx, [yx], [2.0], vol, 1)
        dx, electric_field = charge.calculate_electric_field()
        expected = np.array([-12.95978256, -10.07983088,  -7.1998792,
                             -4.31992752,  -1.43997584, 1.43997584,
                             4.31992752,   7.1998792,   10.07983088,
                             12.95978256])
        assert_almost_equal(electric_field, expected)

    def test_calculate_electrostatic_potential(self):
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        xx, yx, vol = test_density.one_dimensional_density(direction="x")
        charge = analysis.OneDimensionalChargeDensity(xx, [yx], [2.0], vol, 1)
        electrostatic_pot = charge.calculate_electrostatic_potential()[1]
        expected = np.array([-0.0000000e+00,  1.1519807e+01,  2.0159662e+01,
                             2.5919565e+01, 2.8799517e+01,  2.8799517e+01,
                             2.5919565e+01,  2.0159662e+01,
                             1.1519807e+01, -1.7763568e-15])
        assert_almost_equal(electrostatic_pot, expected, decimal=4)
