import sys 

sys.path.append('/home/adam/data/adam/progs/PolyPy/PolyPy/polypy/')

import Generic as ge
import numpy as np
import numpy.testing as npt
from numpy.testing import assert_almost_equal, assert_equal

def test_bin_choose():
    a = ge.bin_choose(4, 2)
    expected = 1
    assert a == expected
    
def test_pbc():
    a, b = ge.pbc(10, 9, 20)
    c, d = ge.pbc(1, 40, 20)
    expected_a = False
    expected_b = 10
    expected_c = True
    expected_d = 41
    assert a == expected_a
    assert b == expected_b
    assert c == expected_c
    assert d == expected_d