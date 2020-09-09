import numpy as np
import os
from polypy import read
from polypy import msd as msd
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'HISTORY')
test_config = os.path.join(os.path.dirname(__file__), 'CONFIG')

class testMSDContainer(unittest.TestCase):

    def test_smooth_msd_data(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        x_data = np.array([1, 2, 3, 4, 5, 1, 2, 3, 4, 5])
        y_data = np.array([5, 4, 3, 2, 1, 5, 4, 3, 2, 1])
        x, y = msds.smooth_msd_data(x_data, y_data)
        predicted_x = np.array([1, 2, 3, 4, 5])
        predicted_y = np.array([5, 4, 3, 2, 1])
        assert_almost_equal(x, predicted_x)
        assert_almost_equal(y, predicted_y)

    def test_xyz_diffusion_coefficient(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        assert_almost_equal(msds.xyz_diffusion_coefficient(), 299.9999999)

    def test_xy_diffusion_coefficient(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        assert_almost_equal(msds.xy_diffusion_coefficient(), 300.0)

    def test_xz_diffusion_coefficient(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        assert_almost_equal(msds.xz_diffusion_coefficient(), 300.0)

    def test_yz_diffusion_coefficient(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        assert_almost_equal(msds.yz_diffusion_coefficient(), 300.0)

    def test_x_diffusion_coefficient(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        assert_almost_equal(msds.x_diffusion_coefficient(), 300.0)

    def test_y_diffusion_coefficient(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        assert_almost_equal(msds.y_diffusion_coefficient(), 300.0)

    def test_z_diffusion_coefficient(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        assert_almost_equal(msds.z_diffusion_coefficient(), 300.0)


class TestMSD(unittest.TestCase):

    def msd_fail_1(self):
        data = read.Config(test_history, ['CA'])
        data.trajectory.timesteps = 1
        with self.assertRaises(ValueError): msd.MSD(data.trajectory)

    def msd_fail_2(self):
        data = read.Config(test_history, ['CA'])
        data.trajectory.atom_name.append('F')
        with self.assertRaises(ValueError): msd.MSD(data.trajectory)

    def test_msd(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        msds = msd_data.msd()
        expected_msd = np.array([3, 12, 27, 48, 75])
        assert_almost_equal(msds.msd, expected_msd)

    def test_calculate_distances(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        trajectories = np.split(data.trajectory.fractional_trajectory, data.trajectory.timesteps)
        x, y = msd_data.calculate_distances(trajectories, 1)
        x = np.asarray(x)**2
        msd_ = np.sum(x, axis=1)
        expected_msd = np.array([3, 12, 27, 3, 12, 27, 48, 75, 108])
        assert_almost_equal(msd_, expected_msd)

    def test_squared_displacements(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        trajectories = np.split(data.trajectory.fractional_trajectory, data.trajectory.timesteps)
        x, y = msd_data.calculate_distances(trajectories, 1)
        x = np.asarray(x)
        msd_data.squared_displacements(x, 1)
        expected_msd = np.array([3, 12, 27, 3, 12, 27, 48, 75, 108])
        assert_almost_equal(msd_data.msd_information.msd, expected_msd)

    def test_three_dimension_square_distance(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        trajectories = np.split(data.trajectory.fractional_trajectory, data.trajectory.timesteps)
        x, y = msd_data.calculate_distances(trajectories, 1)
        x = np.asarray(x) ** 2
        msd_data.three_dimension_square_distance(x, 1)
        expected_msd = np.array([3, 12, 27, 3, 12, 27, 48, 75, 108])
        assert_almost_equal(msd_data.msd_information.msd, expected_msd)

    def test_two_dimension_square_distance(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        trajectories = np.split(data.trajectory.fractional_trajectory, data.trajectory.timesteps)
        x, y = msd_data.calculate_distances(trajectories, 1)
        x = np.asarray(x) ** 2
        msd_data.two_dimension_square_distance(x, 1)
        expected_msd = np.array([2,  8, 18,  2,  8, 18, 32, 50, 72])
        assert_almost_equal(msd_data.msd_information.xymsd, expected_msd)

    def test_one_dimension_square_distance(self):
        data = read.History(test_history, ['CA'])
        msd_data = msd.MSD(data.trajectory)
        trajectories = np.split(data.trajectory.fractional_trajectory, data.trajectory.timesteps)
        x, y = msd_data.calculate_distances(trajectories, 1)
        x = np.asarray(x) ** 2
        msd_data.one_dimension_square_distance(x, 1)
        expected_msd = np.array([ 1,  4,  9,  1,  4,  9, 16, 25, 36])
        assert_almost_equal(msd_data.msd_information.xmsd, expected_msd)

class TestRegionalMSD(unittest.TestCase):

    def test_analyse_trajectory(self):
        pass

    def test_initialise_new_trajectory(self):
        pass

    def test_check_trajectory(self):
        pass

    