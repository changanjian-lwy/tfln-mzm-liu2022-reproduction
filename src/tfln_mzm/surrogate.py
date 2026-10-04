"""Gaussian-process surrogate and constrained acquisition for B05 (numpy/scipy only).

Zero-mean GP on standardized outputs, anisotropic Matern-5/2 kernel, fixed nugget,
hyperparameters by maximum log marginal likelihood with analytic gradients
(Rasmussen & Williams 2006, ch. 4-5). Acquisition is improvement times probability
of feasibility (Gardner et al. 2014) with a known objective. Not a physics module:
it only sees (x, y) pairs produced by the model, never paper data.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve,solve_triangular
from scipy.optimize import minimize
from scipy.stats import norm

S5=np.sqrt(5.0)


def matern52(a,b,ls,var):
    d=(np.asarray(a,float)[:,None,:]-np.asarray(b,float)[None,:,:])/ls
    r=np.sqrt((d*d).sum(-1))
    return var*(1+S5*r+5*r*r/3)*np.exp(-S5*r)


class GP:
    def __init__(self,nugget=1e-6,ls_bounds=(.01,3.),var_bounds=(.05,20.),starts=((.1,1.),(.3,1.),(1.,1.))):
        self.nugget=nugget;self.ls_bounds=ls_bounds;self.var_bounds=var_bounds;self.starts=starts

    def _nll(self,theta,x,z):
        """Negative log marginal likelihood and its gradient in (log ls_1..d, log var)."""
        ls=np.exp(theta[:-1]);var=np.exp(theta[-1]);n=len(z)
        d2=((x[:,None,:]-x[None,:,:])/ls)**2;r=np.sqrt(d2.sum(-1));e=np.exp(-S5*r)
        k=var*(1+S5*r+5*r*r/3)*e
        try:
            c=cho_factor(k+self.nugget*np.eye(n),lower=True)
        except np.linalg.LinAlgError:
            return 1e10,np.zeros_like(theta)
        a=cho_solve(c,z)
        nll=.5*z@a+np.log(np.diag(c[0])).sum()+.5*n*np.log(2*np.pi)
        w=cho_solve(c,np.eye(n))-np.outer(a,a)
        g=[.5*(w*(var*5/3*(1+S5*r)*e*d2[...,j])).sum() for j in range(x.shape[1])]
        g.append(.5*(w*k).sum())
        return float(nll),np.array(g)

    def fit(self,x,y):
        x=np.asarray(x,float);y=np.asarray(y,float)
        self.mu=float(y.mean());self.sd=float(y.std()) or 1.0
        z=(y-self.mu)/self.sd;dim=x.shape[1]
        bounds=[tuple(np.log(self.ls_bounds))]*dim+[tuple(np.log(self.var_bounds))]
        best=None
        for ls0,var0 in self.starts:
            r=minimize(self._nll,np.log([ls0]*dim+[var0]),args=(x,z),jac=True,method='L-BFGS-B',bounds=bounds)
            if best is None or r.fun<best.fun:best=r
        self.ls=np.exp(best.x[:-1]);self.var=float(np.exp(best.x[-1]));self.nll=float(best.fun)
        self.x=x;self.chol=cho_factor(matern52(x,x,self.ls,self.var)+self.nugget*np.eye(len(x)),lower=True)
        self.alpha=cho_solve(self.chol,z)
        return self

    def predict(self,xq):
        """Posterior mean and standard deviation of the latent function, in output units."""
        ks=matern52(xq,self.x,self.ls,self.var)
        v=solve_triangular(self.chol[0],ks.T,lower=True)
        s2=np.maximum(self.var-(v*v).sum(0),1e-12*self.var)
        return ks@self.alpha*self.sd+self.mu,np.sqrt(s2)*self.sd


def log_pof(mu,sd,lower=None,upper=None):
    """log P(y>=lower) or log P(y<=upper) under N(mu, sd^2)."""
    if (lower is None)==(upper is None):
        raise ValueError('Give exactly one of lower or upper')
    return norm.logcdf((mu-lower)/sd) if lower is not None else norm.logcdf((upper-mu)/sd)


def choose(j,lp,evaluated,j_inc):
    """Next index: max log(J-J_inc)+log PoF over unevaluated J>J_inc; max log PoF if no incumbent.

    Ties go to the lowest index. Returns None when no candidate is left.
    """
    j=np.asarray(j,float);lp=np.asarray(lp,float);cand=~np.asarray(evaluated,bool)
    if j_inc is not None:cand&=j>j_inc
    if not cand.any():return None
    score=np.full(len(j),-np.inf)
    score[cand]=lp[cand] if j_inc is None else np.log(j[cand]-j_inc)+lp[cand]
    return int(np.argmax(score))
