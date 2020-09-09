import numpy as np
import os
from polypy import read
from polypy import msd as msd
from polypy import plotting
from polypy.density import Density

import unittest
from numpy.testing import assert_almost_equal


test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')

class TestDensity(unittest.TestCase):

    def test_line_plot(self):   
        data = read.History(test_history, ['CA'])
        test_density = Density(data.trajectory, histogram_size=1.0)
        xx, yx, vol = test_density.one_dimensional_density(direction="x")
        ax = plotting.one_dimensional_density_plot([xx], [ yx], ["Ca"])
        assert_almost_equal(xx, ax.lines[0].get_xydata().T[0] )
        assert_almost_equal(yx, ax.lines[0].get_xydata().T[1] )

    def test_volume_plot(self):  
        x = np.array([1, 2, 3, 4, 5, 6])
        y = np.array([6, 5, 4, 3, 2, 1])
        ax = plotting.volume_plot(x, y)
        assert_almost_equal(x, ax.lines[0].get_xydata().T[0] )
        assert_almost_equal(y, ax.lines[0].get_xydata().T[1] )

    def test_electric_field_plot(self):  
        x = np.array([1, 2, 3, 4, 5, 6])
        y = np.array([6, 5, 4, 3, 2, 1])
        ax = plotting.electric_field_plot(x, y)
        assert_almost_equal(x, ax.lines[0].get_xydata().T[0] )
        assert_almost_equal(y, ax.lines[0].get_xydata().T[1] )

    def test_electrostatic_potential_plot(self):  
        x = np.array([1, 2, 3, 4, 5, 6])
        y = np.array([6, 5, 4, 3, 2, 1])
        ax = plotting.electrostatic_potential_plot(x, y)
        assert_almost_equal(x, ax.lines[0].get_xydata().T[0] )
        assert_almost_equal(y, ax.lines[0].get_xydata().T[1] )

    def test_one_dimensional_charge_density_plot(self):  
        x = np.array([1, 2, 3, 4, 5, 6])
        y = np.array([6, 5, 4, 3, 2, 1])
        ax = plotting.one_dimensional_charge_density_plot(x, y)
        assert_almost_equal(x, ax.lines[0].get_xydata().T[0] )
        assert_almost_equal(y, ax.lines[0].get_xydata().T[1] )

    def test_one_dimensional_density_plot(self):  
        x = np.array([1, 2, 3, 4, 5, 6])
        y = np.array([6, 5, 4, 3, 2, 1])
        ax = plotting.one_dimensional_density_plot([x], [y], ["Ca"])
        assert_almost_equal(x, ax.lines[0].get_xydata().T[0] )
        assert_almost_equal(y, ax.lines[0].get_xydata().T[1] )
    