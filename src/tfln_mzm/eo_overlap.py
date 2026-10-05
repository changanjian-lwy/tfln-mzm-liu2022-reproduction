"""Quasi-static electro-optic overlap on the optical cross-section (SC-05, Track O, O02). Not imported by L0/L1.

The potential is solved on its own enlarged mesh (optical_fem.build_mesh with a larger extent) and E_x is carried to
the optical mesh by point evaluation: grad(phi) of a P2 potential is exactly piecewise linear, so it is represented
without loss on a DG-P1 basis and probed at the optical quadrature points (O02 BOUNDARY section 5). O02 drive: half domain x >= 0 with phi = 0 on x = 0 (rails at +-V/2, gap voltage V), the right rail an
equipotential phi = V/2, all other outer edges natural (zero normal D). eps_DC,LN = diag(eps33, eps11) in (x, y),
because x is crystal Z and y crystal X (optical_fem). P2 potential, E = -grad(phi). Lengths um, fields V/um,
r33 in um/V. Low-frequency constants: data/external/dc_materials/SOURCE.md.

Convention A (ref [10] Eq.(1)(2), Lecture 4 p.22, static.py): Gamma = (g/V) int_LN E_x |E_x,opt|^2 / int |E_x,opt|^2,
VpiL = lam g / (2 n_e^3 r33 Gamma). Convention B: add deps_xx = -n_e^4 r33 <E_x> (element mean) to the LN optical
permittivity and re-solve TE0 at +-V; push-pull arms move by +-dn/dV V, so VpiL = lam / (4 |dn/dV|) (O02 BOUNDARY
section 5 corrects the boundary's lam / (2 |dn/dV|)). Both are per rail section (duty cycle not included).
"""
import numpy as np
from tfln_mzm.optical_fem import LN_NAMES

DC_LN={'clamped':{'eps11':44.0,'eps33':27.9,'r33':31e-6},'unclamped':{'eps11':85.0,'eps33':28.7,'r33':33e-6}}
DC_SUBSTRATE={'quartz':4.5,'fused':3.8}


def dc_permittivity(sec,state,eps_sio2):
    """{region: (eps_xx, eps_yy)} at low frequency; rails are conductors (value unused)."""
    ln=DC_LN[state]
    eps={k:(ln['eps33'],ln['eps11']) for k in LN_NAMES}
    eps.update(clad=(eps_sio2,)*2,box=(eps_sio2,)*2,substrate=(DC_SUBSTRATE[sec.substrate],)*2,air=(1.0,1.0),
               rail_l=(1.0,1.0),rail_r=(1.0,1.0))
    return eps


def solve_potential(mesh,eps_by_region,fixed):
    """Anisotropic Laplace problem div(diag(exx, eyy) grad phi) = 0 with P2 elements.

    fixed: list of (dof indices, value) Dirichlet sets; every other boundary is natural. Returns (basis, phi).
    """
    from skfem import Basis,BilinearForm,ElementTriP0,ElementTriP2,condense,solve
    basis=Basis(mesh,ElementTriP2())
    b0=basis.with_element(ElementTriP0())
    exx,eyy=b0.zeros(),b0.zeros()
    for name,(ex,ey) in eps_by_region.items():
        if name not in mesh.subdomains:continue
        d=b0.get_dofs(elements=name).all();exx[d]=ex;eyy[d]=ey
    if np.any(exx==0) or np.any(eyy==0):raise ValueError('element without a permittivity')

    @BilinearForm
    def a(u,v,w):return w.exx*u.grad[0]*v.grad[0]+w.eyy*u.grad[1]*v.grad[1]
    A=a.assemble(basis,exx=b0.interpolate(exx),eyy=b0.interpolate(eyy))
    x=basis.zeros();D=[]
    for dofs,val in fixed:
        dofs=np.asarray(dofs);x[dofs]=val;D.append(dofs)
    return basis,solve(*condense(A,x=x,D=np.unique(np.concatenate(D))))


def drive_dofs(basis,V,rail='rail_r'):
    """O02 Dirichlet sets: phi = V/2 on every dof of the rail, phi = 0 on the x = 0 edge."""
    f0=basis.mesh.facets_satisfying(lambda p:np.abs(p[0])<1e-9,boundaries_only=True)
    return [(basis.get_dofs(elements=rail).all(),V/2),(basis.get_dofs(facets=f0).all(),0.0)]


def ln_mask(basis):
    """P0 indicator of the LN regions on a basis sharing the mode's mesh."""
    from skfem import ElementTriP0
    b0=basis.with_element(ElementTriP0());mask=b0.zeros()
    for n in LN_NAMES:
        if n in basis.mesh.subdomains:mask[b0.get_dofs(elements=n).all()]=1
    return b0,mask


class FieldX:
    """E_x = -d phi/dx of a P2 potential, evaluable at arbitrary points (exact DG-P1 representation)."""
    def __init__(s,basis,phi):
        from skfem import ElementDG,ElementTriP1
        s.b=basis.with_element(ElementDG(ElementTriP1()))
        s.v=s.b.project(-basis.interpolate(phi).grad[0])
    def at(s,pts):
        """pts: (2, ...) global coordinates -> E_x with shape pts.shape[1:]."""
        flat=pts.reshape(2,-1)
        return (s.b.probes(flat)@s.v).reshape(pts.shape[1:])


def quad_basis(mode,intorder=6):
    """Fresh basis with the mode's element (same dof numbering) and a 6th-order rule (ERRATA E-03)."""
    from skfem import Basis
    return Basis(mode.basis.mesh,mode.basis.elem,intorder=intorder)


def gamma_a(mode,field,gap,V,intorder=6):
    """Convention A overlap; field is a FieldX (possibly on another mesh).

    The integrand E_x,rf |E_x,opt|^2 has degree 5, which the 3-point rule femwell's mode basis inherits would
    under-integrate, so it is integrated on quad_basis (O02 BOUNDARY section 5).
    """
    from skfem import Functional
    bm=quad_basis(mode,intorder)
    ex_rf=field.at(bm.mapping.F(bm.X))
    (et,et_b),_=bm.split(mode.E)
    eo2=np.abs(np.asarray(et_b.interpolate(et))[0])**2
    b0,mask=ln_mask(bm)

    @Functional
    def num(w):return w.mask*w.erf*w.eo2
    @Functional
    def den(w):return w.eo2
    n=num.assemble(bm,mask=b0.interpolate(mask),erf=ex_rf,eo2=eo2)
    return abs(gap/V*n/den.assemble(bm,eo2=eo2))


def element_ex(mesh,field,intorder=6):
    """Element-mean E_x on `mesh` (the optical mesh) from a FieldX, indexed by element."""
    from skfem import Basis,ElementTriP0
    b0=Basis(mesh,ElementTriP0(),intorder=intorder)
    vals=b0.project(field.at(b0.mapping.F(b0.X)))
    return vals[b0.interior_dofs[0]]   # interior_dofs[0][k] is the dof of element k


def deps_xx(mesh,ex_elem,ne,r33,V):
    """Convention B perturbation of eps_xx: -ne^4 r33 E_x in LN elements, zero elsewhere."""
    d=np.zeros(mesh.nelements)
    for n in LN_NAMES:
        if n in mesh.subdomains:
            e=mesh.subdomains[n];d[e]=-ne**4*r33*ex_elem[e]*V
    return d


def vpil_a(lam,gap,ne,r33,gamma):
    """Convention A VpiL per rail section, V*cm (lam, gap um; r33 um/V)."""
    return lam*gap/(2*ne**3*r33*gamma)*1e-4


def vpil_b(lam,dn_dv):
    """Convention B VpiL per rail section, V*cm: push-pull phase difference (4 pi / lam) |dn/dV| V L = pi."""
    return lam/(4*abs(dn_dv))*1e-4
