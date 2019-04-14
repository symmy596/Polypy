import numpy as np
import numpy.testing as npt
import os
from polypy import read_utils as rd_ut
import unittest
from numpy.testing import assert_almost_equal

test_history = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/HISTORY1')
test_config1 = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/CONFIG1')
test_config2 = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/CONFIG2')
test_config3 = os.path.join(os.path.dirname(__file__), 'dlppy_test_data/sample_configs/CONFIG3')
    

test_data_dir='dlppy_test_data/'

class TestReadUtils(unittest.TestCase):
    
    def test_open_file(self):
        file_name = test_config1
        fp_config = rd_ut.open_file(file_name)
        expect = file_name
        assert fp_config[1] == expect