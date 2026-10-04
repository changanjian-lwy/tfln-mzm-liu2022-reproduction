"""A12 module gate: complex characteristic impedance (branch C-01b); independent oracles."""
import unittest
import numpy as np
from scipy.integrate import quad
from tfln_mzm.microwave import propagation, C0
from tfln_mzm.termination import reflection_coefficients, s11
from tfln_mzm.interaction import average_voltage

L=.006;ZG=50.0


def close(a,b,atol=1e-10,rtol=1e-8):
    np.testing.assert_allclose(a,b,atol=atol,rtol=rtol)


def boundary_value_vavg(f,z0,zl,nm,ng,alpha):
    """Solve source/load boundary conditions directly, then integrate V(z)exp(+j*beta_o*z)."""
    g=complex(propagation(f,nm,alpha));bo=2*np.pi*f*ng/C0
    load_row=[0.0,0.0]
    if np.isinf(zl):load_row=[np.exp(-g*L)/z0,-np.exp(g*L)/z0]          # I(L)=0
    else:load_row=[np.exp(-g*L)*(1-zl/z0),np.exp(g*L)*(1+zl/z0)]      # V(L)=ZL I(L)
    vp,vm=np.linalg.solve(np.array([[1+ZG/z0,1-ZG/z0],load_row]),[1.0,0.0])  # Vg=V(0)+Zg I(0), Vg=1
    v=lambda z:(vp*np.exp(-g*z)+vm*np.exp(g*z))*np.exp(1j*bo*z)
    return (quad(lambda z:v(z).real,0,L,epsabs=1e-13,limit=200)[0]+1j*quad(lambda z:v(z).imag,0,L,epsabs=1e-13,limit=200)[0])/L


def physical_z0(zr,alpha,f,nm):
    """G=0 line: Z0=gamma/(j*omega*C) => Im/Re = -alpha/beta."""
    beta=2*np.pi*f*nm/C0
    return zr*(1-1j*alpha/beta)


class TestComplexZ0(unittest.TestCase):
    def test_a_zero_imaginary_part_matches_real_path(self):
        f=np.linspace(1e6,200e9,101)
        for load in [20,50,80,np.inf]:
            close(reflection_coefficients(43+0j,ZG,load),reflection_coefficients(43,ZG,load),atol=1e-15)
            g=propagation(f,2.3,40)
            close(s11(g,L,43+0j,load),s11(g,L,43,load),atol=1e-15)
            close(average_voltage(f,length_m=L,z0=43+0j,zg=ZG,zl=load,n_m=2.3,n_g=2.25,alpha_np_m=40),
                  average_voltage(f,length_m=L,z0=43,zg=ZG,zl=load,n_m=2.3,n_g=2.25,alpha_np_m=40),atol=1e-15)

    def test_b_s11_matches_independent_abcd(self):
        f=np.linspace(1e9,200e9,57)
        for zr in [38.,46.,55.]:
            for load in [0,20,40,80,np.inf]:
                g=propagation(f,2.25,60);z0=zr*(1-.08j)
                a=np.cosh(g*L);b=z0*np.sinh(g*L);c=np.sinh(g*L)/z0
                num,den=(a,c) if np.isinf(load) else (a*load+b,c*load+a)
                close(s11(g,L,z0,load),(num-50*den)/(num+50*den))

    def test_c_vavg_matches_boundary_value_solution(self):
        for f in [2e9,37e9,150e9]:
            for load in [20,80,np.inf]:
                for alpha in [0,90]:
                    z0=physical_z0(44.,alpha,f,2.4) if alpha else 44-6j
                    close(average_voltage(f,length_m=L,z0=z0,zg=ZG,zl=load,n_m=2.4,n_g=2.25,alpha_np_m=alpha),
                          boundary_value_vavg(f,z0,load,2.4,2.25,alpha))

    def test_d_passivity_of_the_g0_branch(self):
        f=np.linspace(.5e9,200e9,400)
        for zr in [37.,45.,55.]:
            for alpha in [5.,50.,120.]:
                z0=physical_z0(zr,alpha,f,2.25)
                for load in [0,20,40,50,80,np.inf]:
                    self.assertTrue(np.all(abs(s11(propagation(f,2.25,alpha),L,z0,load))<=1+1e-12))

    def test_e_invalid_complex_impedance_rejected_and_not_truncated(self):
        for bad in [-1+2j,0+5j,np.nan+1j]:
            with self.assertRaises(ValueError):reflection_coefficients(bad,ZG,50)
        r1,r2=reflection_coefficients(45-5j,ZG,80)
        self.assertNotEqual(np.imag(r1),0);self.assertNotEqual(np.imag(r2),0)


if __name__=='__main__':
    unittest.main()
