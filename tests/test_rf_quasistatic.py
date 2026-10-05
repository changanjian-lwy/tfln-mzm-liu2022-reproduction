"""Track Q (SC-06) quasi-static RF cross-section: exact references, solver and Liu-section geometry."""
import unittest
from collections import OrderedDict
import numpy as np
from scipy.special import ellipk
from tfln_mzm import rf_quasistatic as q

try:
    import femwell  # noqa: F401
    HAVE_FEMWELL=True
except ImportError:
    HAVE_FEMWELL=False


class ExactReferences(unittest.TestCase):
    def test_infinite_ground_modulus(self):
        self.assertAlmostEqual(q.k_cpw(80,20),2/3,places=15)
        self.assertAlmostEqual(q.k_cpw(80,20,1e9),2/3,places=9)
        self.assertLess(q.k_cpw(80,20,50),q.k_cpw(80,20,150))   # narrower grounds -> smaller k -> higher Z0

    def test_symmetric_modulus_and_air_impedance(self):
        k=1/np.sqrt(2);s=2*k;g=1-k                            # a/b = k with infinite grounds
        self.assertAlmostEqual(q.c_cpw_exact(s,g)/(4*q.EPS0),1.0,places=12)
        z_air=1/(q.C_LIGHT*q.c_cpw_exact(80,20))
        kk=2/3;eta0=1/(q.EPS0*q.C_LIGHT)                      # textbook 30 pi uses eta0 ~ 120 pi
        self.assertAlmostEqual(z_air,eta0/4*ellipk(1-kk*kk)/ellipk(kk*kk),places=9)
        self.assertAlmostEqual(q.C_LIGHT*q.inductance(q.c_cpw_exact(80,20)),z_air,places=9)

    def test_loaded_weights(self):
        self.assertAlmostEqual(q.loaded_c(1.0,2.0,3.0,duty=0.9,fn=0.05),0.85*2+0.05*3+0.1*1)


@unittest.skipUnless(HAVE_FEMWELL,"femwell not installed (see requirements-fem-lock.txt)")
class Solver(unittest.TestCase):
    def test_anisotropic_parallel_plates(self):
        from shapely.geometry import LineString,box
        sh=OrderedDict(a=LineString([(0.2,0),(0.2,1)]),b=LineString([(0.8,0),(0.8,1)]),d=box(0,0,1,1))
        mesh,_=q.build_mesh(sh,[],1.0,1.0,far=0.1)
        c,_,_=q.solve_capacitance(mesh,{'d':(3.0,7.0)},[],[('a',1.0),('b',0.0)])
        self.assertAlmostEqual(c/(q.EPS0*3.0/0.6),1.0,places=9)   # exx only; natural top/bottom/outer edges

    def test_cpw_on_anisotropic_half_space(self):
        s,g,wg=10.0,6.0,30.0;X=9.5*(s/2+g+wg)
        mesh,own=q.build_mesh(q.half_space_shapes(s,g,wg,X),['sig','gnd'],2.0,X)
        lines=[('sig',1.0)]+[(k,1.0) for k in own['sig']]+[('gnd',0.0)]+[(k,0.0) for k in own['gnd']]
        c,_,_=q.solve_capacitance(mesh,{'above':(1.0,1.0),'below':(27.9,44.0)},[],lines)
        self.assertLess(abs(2*c/q.c_cpw_exact(s,g,wg,np.sqrt(27.9*44.0))-1),2e-3)


class LiuSection(unittest.TestCase):
    def test_kinds_and_windows(self):
        X=2000.0
        for etch in ('full','gap'):
            for kind in ('U','H','N'):
                for win in (None,4.5):
                    for t in (0.0,1.5):
                        sec=q.RFSection(t_bcb=t,kind=kind,etch=etch,window=win)
                        sh,cond,thin=q.liu_shapes(sec,X)
                        self.assertEqual(cond[:2],['sig','gnd'])
                        self.assertEqual(len(cond),2 if kind=='U' else 4)
                        if t>0:
                            self.assertTrue(sh['bcb'].is_valid)
                            inside=sh['bcb'].intersection(q.box(sec.xa-sec.W+1e-6,0,sec.xa+sec.W-1e-6,X))
                            self.assertLess(inside.area,1e-9)                       # no BCB in the window
                        for name in cond:self.assertTrue(sh[name].is_valid)
                        if kind!='U':
                            gap=sh['head_g'].bounds[0]-sh['head_s'].bounds[2]
                            self.assertAlmostEqual(gap,q.HEAD_GAP,places=9)
                        y_el=sh['sig'].bounds[1]
                        base=(0.3 if etch=='full' else 0.6)+0.1
                        self.assertAlmostEqual(y_el,base+t,places=9)

    def test_neck_reaches_electrode(self):
        for etch in ('full','gap'):
            for win in (None,4.5):
                sec=q.RFSection(kind='N',etch=etch,window=win)
                sh,_,_=q.liu_shapes(sec,2000.0)
                self.assertLess(sh['head_s'].distance(sh['sig']),1e-9)
                self.assertLess(sh['head_g'].distance(sh['gnd']),1e-9)


if __name__=='__main__':
    unittest.main()
