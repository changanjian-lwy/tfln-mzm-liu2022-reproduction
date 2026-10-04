import unittest
import numpy as np
from scipy.stats import norm
from tfln_mzm.surrogate import GP,choose,log_pof,matern52


def smooth(x):
    return np.sin(3*x[:,0])+np.cos(2*x[:,1])+x[:,0]*x[:,1]


class TestSurrogate(unittest.TestCase):
    def setUp(self):
        self.x=np.random.default_rng(0).random((30,2))

    def test_kernel_value_and_symmetry(self):
        k=matern52(self.x[:5],self.x[:5],np.array([.3,.7]),2.0)
        self.assertTrue(np.allclose(k,k.T));self.assertTrue(np.allclose(np.diag(k),2.0))
        r=np.sqrt(5)*0.5
        self.assertAlmostEqual(matern52([[0,0]],[[.15,0]],np.array([.3,1.]),1.0)[0,0],(1+r+r*r/3)*np.exp(-r))

    def test_nll_gradient_matches_finite_difference(self):
        gp=GP();y=smooth(self.x);z=(y-y.mean())/y.std();t=np.log([.4,.6,1.3])
        _,g=gp._nll(t,self.x,z)
        for i in range(3):
            h=np.zeros(3);h[i]=1e-6
            fd=(gp._nll(t+h,self.x,z)[0]-gp._nll(t-h,self.x,z)[0])/2e-6
            self.assertLess(abs(fd-g[i]),1e-5*max(1,abs(fd)))

    def test_interpolates_training_data(self):
        y=smooth(self.x);m,s=GP().fit(self.x,y).predict(self.x)
        self.assertLess(np.abs(m-y).max(),1e-3*y.std());self.assertLess(s.max(),1e-2*y.std())

    def test_recovers_smooth_function(self):
        g=np.stack(np.meshgrid(np.linspace(0,1,21),np.linspace(0,1,21)),-1).reshape(-1,2)
        m,s=GP().fit(self.x,smooth(self.x)).predict(g);y=smooth(g)
        self.assertLess(np.sqrt(np.mean((m-y)**2)),.05*y.std())
        self.assertGreater(np.mean(np.abs(m-y)<=3*s),.9)

    def test_log_pof(self):
        self.assertAlmostEqual(float(log_pof(0.,1.,lower=0.)),np.log(.5))
        self.assertAlmostEqual(float(log_pof(1.,2.,upper=3.)),float(norm.logcdf(1.)))
        with self.assertRaises(ValueError):log_pof(0.,1.)

    def test_choose_rules(self):
        j=np.array([1.,3.,3.,2.,5.]);lp=np.array([0.,-1.,-1.,0.,-50.]);ev=np.array([True,False,False,False,False])
        self.assertEqual(choose(j,lp,ev,None),3)
        self.assertEqual(choose(j,lp,ev,2.0),1)
        self.assertEqual(choose(j,np.array([0.,-1.,-1.,0.,-1.]),ev,2.0),4)
        self.assertIsNone(choose(j,lp,np.ones(5,bool),None))
        self.assertIsNone(choose(j,lp,ev,5.0))


if __name__=='__main__':
    unittest.main()
