import TrajectoryAnalysis as ta
import numpy as np
import numpy.testing as npt
from numpy.testing import assert_almost_equal, assert_equal

def test_distances():
    distance = ta.distances(10, 5)
    expect = 5
    assert distance == expect

def test_msd_stats():
    data = np.genfromtxt("test_data/MSD.txt", dtype="float")
    a, b, c, d = ta.msd_stats(data[:,0], data[:,1], data[:,2], data[:,3], data[:,4])
    assert_almost_equal(a, 100.0)
    assert_almost_equal(b, 10.0)
    assert_almost_equal(c, 10.0)
    assert_almost_equal(d, 10.0)

def test_diffusion_coefficient():
    a, b, c, d = ta.diffusion_coefficient(24, 25, 14, 73)
    assert_almost_equal(a, 40)
    assert_almost_equal(b, 125)
    assert_almost_equal(c, 70)
    assert_almost_equal(d, 365)
