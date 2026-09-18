"""Independent physics limits B01–B09; no digitized-curve fit."""
import unittest
import numpy as np
from scipy.integrate import quad
from tfln_mzm.microwave import propagation, db_per_mm_to_np_per_m, C0
from tfln_mzm.termination import reflection_coefficients, s11
from tfln_mzm.interaction import average_voltage
from tfln_mzm.response import eo_response_db, magnitude_db

P=dict(length_m=.006,z0=50,zg=50,zl=50,n_m=2.25,n_g=2.25,alpha_np_m=0,vg=1)
F=np.linspace(0,200e9,301)

def av(f=F, **kw):
    return average_voltage(f,**(P|kw))

def close(a,b):
    np.testing.assert_allclose(a,b,atol=1e-10,rtol=1e-8)

class TestTravelingWave(unittest.TestCase):
    def test_B01_dc_divider(self):
        for z in [35,50,65]:
            for load in [20,40,50,80]:
                close(av(0,z0=z,zl=load),load/(50+load))

    def test_B02_perfect_matching(self):
        close(av(),.5)
        close(eo_response_db(av(),av(1e6)),0)
        close(s11(propagation(F,2.25,0),.006,50,50),0)

    def test_B03_matched_loss(self):
        for alpha in [20,100]:
            close(abs(av(alpha_np_m=alpha)),.5*(-np.expm1(-alpha*.006))/(alpha*.006))

    def test_B04_velocity_mismatch(self):
        close(abs(av(n_m=2.4)),.5*abs(np.sinc(F*.15*.006/C0)))

    def test_B05_open_short(self):
        for load,r,voltage in [(np.inf,1,1),(0,-1,0)]:
            close(reflection_coefficients(50,50,load)[1],r)
            close(av(0,zl=load),voltage)
            close(abs(s11(propagation(F,2.25,0),.006,50,load)),1)

    def test_B06_zero_reference_rejected(self):
        with self.assertRaises(ValueError):
            eo_response_db(av(zl=0),av(0,zl=0))

    def test_B07_short_length(self):
        for load in [20,50,80,np.inf]:
            expected=1 if np.isinf(load) else load/(50+load)
            close(av(length_m=1e-12,z0=42,zl=load,alpha_np_m=20),expected)

    def test_B08_passivity_and_independent_abcd(self):
        for load in [0,20,50,80,np.inf]:
            g=propagation(F,2.35,20)
            actual=s11(g,.006,43,load)
            self.assertTrue(np.all(abs(actual)<=1+1e-12))
            a=np.cosh(g*.006);b=43*np.sinh(g*.006);c=np.sinh(g*.006)/43
            # Direct ABCD terminated input reflection; avoids infinite Zin.
            num,den=(a,c) if np.isinf(load) else (a*load+b,c*load+a)
            close(actual,(num-50*den)/(num+50*den))
        self.assertEqual(float(magnitude_db(0)),-np.inf)

    def test_B09_printed_equation_quadrature(self):
        for f in [0,1e6,37e9,200e9]:
            for load in [20,80,np.inf]:
                for alpha in [0,100]:
                    L=.006;z0=43;zg=50
                    be=2*np.pi*f*2.4/C0-1j*alpha
                    bo=2*np.pi*f*2.25/C0
                    # Independent literal Eq.(1), rather than the stable helper.
                    r1=(z0-zg)/(z0+zg)
                    r2=1 if np.isinf(load) else (load-z0)/(load+z0)
                    def integrand(t):
                        x=L*t
                        return (1+r1)*np.exp(1j*bo*L)*(np.exp(1j*(be-bo)*x)+r2*np.exp(-1j*(be+bo)*x))/(2*(np.exp(1j*be*L)+r1*r2*np.exp(-1j*be*L)))
                    expected=quad(lambda t:integrand(t).real,0,1,epsabs=1e-12,epsrel=1e-12)[0]+1j*quad(lambda t:integrand(t).imag,0,1,epsabs=1e-12,epsrel=1e-12)[0]
                    close(av(f,z0=z0,zl=load,n_m=2.4,alpha_np_m=alpha),expected)

    def test_loss_units_and_invalid_inputs(self):
        close(db_per_mm_to_np_per_m(1),1000*np.log(10)/20)
        for kw in [dict(length_m=0),dict(alpha_np_m=-1),dict(z0=0),dict(zl=-1),dict(n_g=0)]:
            with self.assertRaises(ValueError):av(**kw)
        with self.assertRaises(ValueError):av(-1)
