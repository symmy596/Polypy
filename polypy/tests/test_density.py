import numpy as np
import os
from polypy import read as rd
from polypy import density as dens
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')
test_config = os.path.join(os.path.dirname(__file__), 'CONFIG')
test_archive = os.path.join(os.path.dirname(__file__), 'ARCHIVE')


class TestUtils(unittest.TestCase):

    def test_one_dimensional_density(self):
        data = rd.read_history(test_history, ["CA"])
        test_density = dens.Density(data)
        x, y = test_density.one_dimensional_density(Bin=1.00)
        predicted_x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        predicted_y = np.zeros(10) + 1
        assert_almost_equal(x, predicted_x)
        assert_almost_equal(y, predicted_y)