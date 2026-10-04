"""Coplanar-waveguide cross-section FEM on femwell (SC-01, layer L2). Not imported by the L0/L1 code.

femwell (GPL-3.0, pinned in requirements-fem-lock.txt) is used as a dependency; no femwell code is copied.
Method follows femwell's RF CPW tutorial: one quasi-TEM eigenmode of the x>=0 half (the x=0 symmetry
plane and the outer edges are natural, i.e. PMC, boundaries); Z0 from the power-current definition
p0/|i0|^2 and R, L, G, C from field integrals (Marks & Williams 1992). Unlike the tutorial, i0 is the
conduction current integral of sigma*E_z over the signal strip. Mesh lengths are micrometres; returned
quantities are SI per full structure. Below ~50 MHz the eigenproblem is unreliable on the Tuncer CPW.

Mesh scale m multiplies every mesh size (resolution, SizeMax, default maximum) of the tutorial's "fine"
mesh; transition distances stay fixed. Domain scale d multiplies the tutorial's sweep domain
(half-width 40.5 um, 27 um above and below); a ground plane is clipped to the domain.
"""
from collections import OrderedDict
from dataclasses import dataclass
import time
import numpy as np
from scipy.constants import epsilon_0, mu_0, speed_of_light
from shapely.geometry import LineString, box
from shapely.ops import clip_by_rect, linemerge

CONTOUR='sig_contour'


@dataclass(frozen=True)
class CPW:
    """Symmetric CPW on a substrate; metal from y=0 to y=t_metal, substrate below y=0, air above. Lengths in um."""
    w_sig:float
    gap:float
    w_gnd:float
    t_metal:float
    sigma:float      # metal conductivity, S/m
    eps_sub:float    # substrate relative permittivity (lossless)


def tutorial_domain(cpw,d):
    """Half-width, height above, depth below (um) of femwell's sweep domain scaled by d."""
    span=cpw.w_sig+2*cpw.gap
    return 1.5*span*d,span*d,span*d


def half_geometry(cpw,d):
    """Ordered shapes for the x>=0 half; earlier entries override later ones when meshed."""
    W,top,bot=tutorial_domain(cpw,d)
    w2,t=cpw.w_sig/2,cpw.t_metal
    g0=w2+cpw.gap;g1=min(g0+cpw.w_gnd,W)
    sig_full=box(-w2,0,w2,t);gnd=box(g0,0,g1,t)
    off=min(t/20,cpw.w_sig/10)
    ring=LineString(sig_full.buffer(off,join_style='bevel').exterior.coords)
    contour=clip_by_rect(ring,0,-bot,W,top)
    if contour.geom_type=='MultiLineString':contour=linemerge(contour)
    gring=LineString(gnd.buffer(min(t/20,cpw.w_gnd/10),join_style='bevel').exterior.coords)
    gring=clip_by_rect(gring,0,-bot,W,top)
    if gring.geom_type=='MultiLineString':gring=linemerge(gring)
    gap_line=LineString([(0,t/2),((g0+g1)/2,t/2)])
    return OrderedDict([(CONTOUR,contour),('gnd_contour',gring),('gap_line',gap_line),
        ('metal_sig',box(0,0,w2,t)),('metal_gnd',gnd),('air',box(0,0,W,top)),('substrate',box(0,-bot,W,0))]),off


def resolutions(m):
    """femwell tutorial 'fine' sizes scaled by m (gnd_contour carries no size, as in the tutorial)."""
    return {CONTOUR:{'resolution':.1*m,'distance':10},'gap_line':{'resolution':.2*m,'distance':10},
            'metal_sig':{'resolution':.1*m,'distance':.1,'SizeMax':.1*m},
            'metal_gnd':{'resolution':.5*m,'distance':.2,'SizeMax':.5*m},
            'air':{'resolution':10*m,'distance':.1},'substrate':{'resolution':10*m,'distance':1}}


def build_mesh(cpw,m,d):
    from femwell.mesh import mesh_from_OrderedDict
    from skfem.io.meshio import from_meshio
    shapes,_=half_geometry(cpw,d)
    return from_meshio(mesh_from_OrderedDict(shapes,resolutions(m),default_resolution_max=100*m))


def solve(cpw,mesh,f,order=1,n_guess=None):
    """Quasi-TEM mode at frequency f (Hz). Returns SI quantities for the full structure.

    order is the femwell element order (1 or 2). n_guess (may be complex) sets the shift-invert target
    k0^2 n_guess^2; None keeps femwell's default shift 1.1 k0^2 max(eps). Defaults reproduce A17.
    """
    from femwell.maxwell.waveguide import compute_modes
    from skfem import Basis,ElementTriP0,Functional
    from skfem.helpers import inner
    omega=2*np.pi*f
    basis0=Basis(mesh,ElementTriP0(),intorder=4)
    eps=basis0.ones(dtype=complex)
    eps[basis0.get_dofs(elements='air')]=1
    eps[basis0.get_dofs(elements='substrate')]=cpw.eps_sub
    for name in ('metal_sig','metal_gnd'):
        eps[basis0.get_dofs(elements=name)]=1-1j*cpw.sigma/(omega*epsilon_0)
    t0=time.perf_counter()
    mode=compute_modes(basis0,eps,wavelength=speed_of_light/f*1e6,num_modes=1,metallic_boundaries=False,
                       order=order,n_guess=n_guess)[0]
    (et,et_b),(ez,ez_b)=mode.basis.split(mode.E)
    (ht,ht_b),(hz,hz_b)=mode.basis.split(mode.H)

    @Functional(dtype=complex)
    def power(w):  # z-component of E_t x conj(H_t)
        return w.e[0]*np.conj(w.h[1])-w.e[1]*np.conj(w.h[0])
    p0=power.assemble(mode.basis,e=et_b.interpolate(et),h=ht_b.interpolate(ht))

    # i0 = conduction current through the signal strip (A17 boundary addendum); the tutorial's contour
    # integral of H converged poorly in the DC-resistance test and is kept only as a diagnostic ratio.
    mask=basis0.zeros();mask[basis0.get_dofs(elements='metal_sig')]=1

    @Functional(dtype=complex)
    def conduction(w):return w.mask*w.ez
    i0=cpw.sigma*1e-6*conduction.assemble(mode.basis,mask=mode.basis_epsilon_r.interpolate(mask),ez=ez_b.interpolate(ez))
    centre=np.array([0.0,cpw.t_metal/2])

    @Functional(dtype=complex)
    def contour(w):  # H_t . tangent, normal oriented away from the strip so facets cannot cancel
        s=np.sign((w.x[0]-centre[0])*w.n[0]+(w.x[1]-centre[1])*w.n[1])
        return s*(-w.n[1]*w.h[0]+w.n[0]*w.h[1])
    fb=ht_b.boundary(facets=mesh.boundaries[CONTOUR])
    i_contour=contour.assemble(fb,h=fb.interpolate(ht))
    v0=p0/np.conj(i0)

    epsF=mode.basis_epsilon_r.interpolate(eps*epsilon_0*1e-6)  # F/um
    muH=mu_0*1e-6                                              # H/um
    fields=dict(eps=epsF,et=et_b.interpolate(et),ez=ez_b.interpolate(ez),ht=ht_b.interpolate(ht),hz=hz_b.interpolate(hz))

    @Functional(dtype=complex)
    def c_form(w):return np.real(w.eps)*inner(w.et,np.conj(w.et))-muH*inner(w.hz,np.conj(w.hz))

    @Functional(dtype=complex)
    def l_form(w):return muH*inner(w.ht,np.conj(w.ht))-np.real(w.eps)*inner(w.ez,np.conj(w.ez))

    @Functional(dtype=complex)
    def g_form(w):return -np.imag(w.eps)*inner(w.et,np.conj(w.et))

    @Functional(dtype=complex)
    def r_form(w):return -np.imag(w.eps)*inner(w.ez,np.conj(w.ez))
    v2,i2=abs(v0)**2,abs(i0)**2
    # half structure -> full: series quantities halve, shunt quantities double; per um -> per m
    C=2*c_form.assemble(mode.basis,**fields).real/v2*1e6
    G=2*omega*g_form.assemble(mode.basis,**fields).real/v2*1e6
    L=l_form.assemble(mode.basis,**fields).real/i2/2*1e6
    R=omega*r_form.assemble(mode.basis,**fields).real/i2/2*1e6
    solve_s=time.perf_counter()-t0
    k=complex(mode.k)*1e6                     # 1/m, field ~ exp(-j k z)
    gamma_fem=1j*k
    zs,ys=R+1j*omega*L,G+1j*omega*C
    neff=complex(mode.n_eff)
    alpha=gamma_fem.real
    return {'f_Hz':f,'neff':neff,'nm':neff.real,'alpha_Np_per_m':alpha,'alpha_dB_per_cm':20/np.log(10)*alpha/100,
            'Z0_pi':complex(p0/i2)/2,'Z0_rlgc':complex(np.sqrt(zs/ys)),'gamma_fem':gamma_fem,'gamma_rlgc':complex(np.sqrt(zs*ys)),
            'R':R,'L':L,'G':G,'C':C,'contour_to_conduction_current':float(abs(i_contour)/abs(i0)),
            'n_elements':int(mesh.nelements),'n_dofs':int(mode.basis.N),'solve_s':solve_s}
