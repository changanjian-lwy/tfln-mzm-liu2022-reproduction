import unittest
import numpy as np
from scipy.optimize import brentq
from tfln_mzm.response import bandwidth_3db,eo_response_db
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import C0

class TestBandwidth(unittest.TestCase):
    def test_B10_sinc_oracle_and_refinement(self):
        f0=1e6;L=.006;dn=.15
        p=dict(length_m=L,z0=50,zg=50,zl=50,n_m=2.4,n_g=2.25,alpha_np_m=0)
        ref=average_voltage(f0,**p)
        oracle=brentq(lambda f:20*np.log10(abs(np.sinc(f*dn*L/C0)/np.sinc(f0*dn*L/C0)))+3,1e9,250e9)
        crossings=[]
        for count in [2001,4001]:
            f=np.linspace(f0,200e9,count)
            r=bandwidth_3db(f,eo_response_db(average_voltage(f,**p),ref))
            crossings.append(r['bandwidth_hz'])
            self.assertLess(abs(r['bandwidth_hz']-oracle),.1e9)
        self.assertLess(abs(crossings[1]-crossings[0]),.1e9)

    def test_B11_censored(self):
        r=bandwidth_3db([1,2,3],[0,-1,-2])
        self.assertIsNone(r['bandwidth_hz']);self.assertTrue(r['censored'])
        self.assertEqual(r['lower_bound_hz'],3)

    def test_B12_first_crossing_not_peak_relative(self):
        r=bandwidth_3db([1,2,3,4,5],[0,2,-4,0,-6])
        self.assertAlmostEqual(r['bandwidth_hz'],2+5/6)
        self.assertEqual(len(r['events']),2)

    def test_contact_and_infinite_null(self):
        self.assertEqual(bandwidth_3db([1,2,3],[0,-3,0])['status'],'threshold_contact')
        r=bandwidth_3db([1,2],[0,-np.inf])
        self.assertFalse(r['censored']);self.assertIsNone(r['bandwidth_hz'])
        self.assertEqual(r['status'],'bracket_requires_refinement')

    def test_bad_grid_and_reference(self):
        for f,y in [([1,1],[0,-4]),([1,2],[1,-4]),([1,2],[0,np.nan]),([1],[0])]:
            with self.assertRaises(ValueError):bandwidth_3db(f,y)
