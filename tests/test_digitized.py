import unittest
import numpy as np
from tfln_mzm.digitized import interpolate_supported

class TestDigitizedSupport(unittest.TestCase):
    def test_no_extrapolation_or_long_gap(self):
        values,mask=interpolate_supported(np.array([-1,0,.5,1,3,5,6]),[0,1,5],[0,2,10],1)
        np.testing.assert_array_equal(mask,[False,True,True,True,False,True,False])
        np.testing.assert_allclose(values[mask],[0,1,2,10])
        self.assertTrue(np.all(np.isnan(values[~mask])))

    def test_invalid_order(self):
        with self.assertRaises(ValueError):interpolate_supported([0],[0,0],[1,2],1)
