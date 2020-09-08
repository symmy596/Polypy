import numpy as np
from polypy import utils as ut
from polypy import read
import unittest
from numpy.testing import assert_almost_equal
import os


test_history = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/HISTORY1')


class TestUtils(unittest.TestCase):

    def test_pbc(self):
        a, b = ut.pbc(0.5, 0.6)
        c, d = ut.pbc(0.1, 0.9)
        expected_a = False
        expected_b = 10
        expected_c = True
        expected_d = 41
        assert a == expected_a
        assert b == 0.5
        assert c == expected_c
        assert d == 1.1
