"""Validate the audit's phase-independent bound with an independent phase sweep."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_a05.py'
spec=importlib.util.spec_from_file_location('audit',p)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

class TestAuditBounds(unittest.TestCase):
    def test_full_phase_sweep_is_bounded(self):
        phase=np.linspace(0,2*np.pi,20001)
        for z in [38.,50.,56.]:
            for load in [20.,40.,50.,80.]:
                for loss in [0.,.3,.8]:
                    low,high=audit.bounds(z,loss,load)
                    r=(z-50)/(z+50)
                    # Independent amplitude attenuation via dB voltage ratio.
                    q=(load-z)/(load+z)*10**(-2*loss*6/20)*np.exp(1j*phase)
                    m=abs((r+q)/(1+r*q))
                    self.assertTrue(np.all(m>=low-1e-12))
                    self.assertTrue(np.all(m<=high+1e-12))
                    self.assertAlmostEqual(float(m.min()),float(low),places=9)
                    self.assertAlmostEqual(float(m.max()),float(high),places=9)

    def test_lossless_reference_shift_preserves_magnitude(self):
        s=np.array([0,.1+.2j,-.5j])
        np.testing.assert_allclose(abs(s*np.exp(-2j*np.array([1,3,8]))),abs(s),atol=1e-14)
