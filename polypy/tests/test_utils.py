import unittest
import numpy as np
from polypy import utils as ut
from numpy.testing import assert_almost_equal


class TestUtils(unittest.TestCase):

    def test_pbc_1(self):
        a, b = ut.pbc(0.5, 0.6)
        c, d = ut.pbc(0.1, 0.9)
        expected_a = False
        expected_c = True
        assert a == expected_a
        assert b == 0.5
        assert c == expected_c
        assert d == 1.1

    def test_pbc_2(self):
        a, b = ut.pbc(5.5, 4.6)
        c, d = ut.pbc(5.1, 5.9)
        expected_a = True
        expected_c = True
        assert a == expected_a
        assert b == 4.5
        assert c == expected_c
        assert d == 6.1

    def test_calculate_rcplvs(self):
        lv = np.array([[3, 0, 0], [0, 5, 0], [0, 0, 2]])
        rclpvs, lengths = ut.calculate_rcplvs(lv)
        expected_rcplvs = np.array([[0.3333333, 0, 0],
                                    [0, 0.2, 0],
                                    [0, 0, 0.5]])
        expected_lengths = np.array([3, 5, 2])
        assert_almost_equal(rclpvs, expected_rcplvs)
        assert_almost_equal(lengths, expected_lengths)

    def test_cart_2_frac_1(self):
        coord = np.array([5, 5, 5])
        lv = np.array([[3, 0, 0], [0, 5, 0], [0, 0, 2]])
        rclpvs, lengths = ut.calculate_rcplvs(lv)
        coords = ut.cart_2_frac(coord, lengths, rclpvs)
        assert_almost_equal(coords, np.array([0.66666667, 0.0, 0.5]))

    def test_cart_2_frac_2(self):
        coord = np.array([[5, 5, 5], [8, 8, 8]])
        lv = np.array([[3, 0, 0], [0, 5, 0], [0, 0, 2]])
        rclpvs, lengths = ut.calculate_rcplvs(lv)
        coords = ut.cart_2_frac(coord, lengths, rclpvs)
        assert_almost_equal(coords, np.array([[0.66666667, 0.0, 0.5],
                                              [0.66666667, 0.6, 0.0]]))
