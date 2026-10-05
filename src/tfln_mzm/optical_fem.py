"""Optical cross-section eigenmodes on femwell (SC-05, Track O). Not imported by the L0/L1 code.

femwell (GPL-3.0, pinned in requirements-fem-lock.txt) is used as a dependency; no femwell code is copied.
Lengths are micrometres. Axes: x horizontal = crystal Z of the X-cut film, y vertical = crystal X,
propagation along crystal Y, so eps_LN = diag(ne^2, no^2, no^2) in femwell's (x, y, z) order.
Materials are read from the refractiveindex.info files in data/external/optical_materials (SOURCE.md).
Lossy media use femwell's sign convention (eps'' < 0, as in cpw_fem), so gold is (n - jk)^2.

Cross-section (O01 BOUNDARY section 1): y = 0 is the bottom of the LN film; bonding SiO2 below it, then
the substrate; LN slab of thickness t_film - etch across the domain; trapezoidal ridge on top; conformal
PECVD SiO2 of thickness t_clad over the LN; air above; optional gold T-rails on the cladding, inner edges
at x = +-gap/2. Mesh scale m multiplies every mesh size; transition distances stay fixed.
"""
from collections import OrderedDict
from dataclasses import dataclass,replace
from functools import lru_cache
from pathlib import Path
import re,time
import numpy as np
from scipy.optimize import brentq
from shapely.geometry import LineString,Polygon,box
from shapely.ops import unary_union

DATA=Path(__file__).resolve().parents[2]/'data/external/optical_materials'
LN_FILES={'CLN':('LiNbO3_Zelmon-e.yml','LiNbO3_Zelmon-o.yml'),'MgO':('MgO-LiNbO3_Zelmon-e.yml','MgO-LiNbO3_Zelmon-o.yml')}
SILICA='SiO2_Malitson.yml'
SUBSTRATE_FILES={'quartz':'SiO2-quartz_Ghosh-o.yml','fused':SILICA}
GOLD='Au_Olmon-ev.yml'
LN_NAMES=('ridge','slab_near','slab')


@lru_cache(maxsize=None)
def _load(name):
    """(kind, coefficients or table, wavelength range in um) from a refractiveindex.info YAML file."""
    text=(DATA/name).read_text()
    kind=re.search(r'type:\s*(.+)',text).group(1).strip()
    if kind.startswith('formula'):
        coef=np.array([float(v) for v in re.search(r'coefficients:\s*(.+)',text).group(1).split()])
        lo,hi=(float(v) for v in re.search(r'wavelength_range:\s*(.+)',text).group(1).split())
        return kind,coef,(lo,hi)
    rows=[r.split() for r in text.split('data: |',1)[1].splitlines()]
    tab=np.array([[float(v) for v in r] for r in rows if len(r)==3])
    return kind,tab,(tab[0,0],tab[-1,0])


def index(name,lam):
    """Complex refractive index n - jk (k >= 0) of a material file at vacuum wavelength lam (um)."""
    kind,c,(lo,hi)=_load(name)
    if not lo<=lam<=hi:raise ValueError(f'{name}: {lam} um outside {lo}-{hi} um')
    if kind=='tabulated nk':return complex(np.interp(lam,c[:,0],c[:,1]),-np.interp(lam,c[:,0],c[:,2]))
    l2=lam*lam;s=1+c[0]
    for i in range(1,len(c),2):
        s+=c[i]*l2/(l2-(c[i+1]**2 if kind=='formula 1' else c[i+1]))
    return complex(np.sqrt(s))


@dataclass(frozen=True)
class Section:
    """One cross-section configuration (O01 BOUNDARY section 1 table). Lengths in um, angle in degrees."""
    theta_deg:float=75.0        # sidewall angle from horizontal
    width:float=1.0
    width_ref:str='top'         # which ridge width equals `width`: 'top' or 'bottom'
    t_film:float=0.6
    etch:float=0.3
    t_clad:float=0.1
    t_box:float=2.0
    gap:float=3.0
    clad_dn:float=0.0           # added to the PECVD cladding index only
    substrate:str='quartz'      # 'quartz' (Ghosh n_o) or 'fused' (Malitson)
    rail_t:float|None=0.5       # T-rail thickness; None = no T-rails
    rail_w:float=2.0
    ln:str='CLN'

    def widths(s):
        """(top, bottom) ridge widths."""
        run=0.0 if s.theta_deg>=90 else s.etch/np.tan(np.radians(s.theta_deg))
        top,bot=(s.width,s.width+2*run) if s.width_ref=='top' else (s.width-2*run,s.width)
        if top<=0:raise ValueError('ridge top width <= 0')
        return top,bot


def domain(sec,d):
    """(x half-width, y bottom, y top) for domain scale d (O01 BOUNDARY: 4d, substrate d, air to 0.6 + 2.1d)."""
    return 4*d,-sec.t_box-d,sec.t_film+2.1*d


def geometry(sec,d=1.0,half=False,extent=None):
    """Ordered shapes; earlier entries override later ones when meshed (femwell mesh_from_OrderedDict).

    half=True keeps x >= 0 only (O01 BOUNDARY section 5): with PEC on every outer edge, x = 0 is a PEC
    symmetry plane, exact for TE0 of this mirror-symmetric section (its E_y and E_z are odd in x).
    extent=(X, ybot, ytop) overrides domain(sec, d) (O02 uses a larger electrostatic domain).
    """
    X,ybot,ytop=extent if extent is not None else domain(sec,d)
    ys=sec.t_film-sec.etch
    top,bot=sec.widths()
    ridge=Polygon([(-bot/2,ys),(bot/2,ys),(top/2,sec.t_film),(-top/2,sec.t_film)])
    slab=box(-X,0,X,ys)
    lnu=unary_union([slab,ridge])
    clad=lnu.buffer(sec.t_clad,join_style='mitre').intersection(box(-X,ys,X,ytop)).difference(lnu)
    shapes=OrderedDict()
    if sec.rail_t is not None:
        g=sec.gap/2;y0=ys+sec.t_clad;y1=y0+sec.rail_t
        shapes['rail_face_l']=LineString([(-g,y0),(-g,y1)]);shapes['rail_face_r']=LineString([(g,y0),(g,y1)])
        shapes['rail_l']=box(-g-sec.rail_w,y0,-g,y1);shapes['rail_r']=box(g,y0,g+sec.rail_w,y1)
    shapes['ridge']=ridge
    shapes['clad']=clad
    shapes['slab_near']=box(-bot/2-1,0,bot/2+1,ys)
    shapes['slab']=slab
    shapes['air']=box(-X,ys,X,ytop)
    shapes['box']=box(-X,-sec.t_box,X,0)
    shapes['substrate']=box(-X,ybot,X,-sec.t_box)
    if half:
        keep=box(0,ybot,X,ytop)
        shapes=OrderedDict((k,v.intersection(keep)) for k,v in shapes.items() if not k.endswith('_l'))
    return shapes


def resolutions(m,rails=True,half=False):
    """O01 BOUNDARY section 1 mesh sizes at scale m (distances fixed)."""
    r={'ridge':{'resolution':.02*m,'distance':.5},'clad':{'resolution':.02*m,'distance':.5},
       'slab_near':{'resolution':.04*m,'distance':1}}
    if rails:
        for side in ('r',) if half else ('l','r'):
            r.update({f'rail_face_{side}':{'resolution':.01*m,'distance':.1},f'rail_{side}':{'resolution':.05*m,'distance':.1}})
    return r


def build_mesh(sec,m,d=1.0,half=False,extent=None):
    from femwell.mesh import mesh_from_OrderedDict
    from skfem.io.meshio import from_meshio
    return from_meshio(mesh_from_OrderedDict(geometry(sec,d,half,extent),resolutions(m,sec.rail_t is not None,half),
                                             default_resolution_max=.4*m))


def permittivity(sec,lam):
    """{region: (eps_xx, eps_yy, eps_zz)} at wavelength lam (um)."""
    fe,fo=LN_FILES[sec.ln]
    ne,no=index(fe,lam).real,index(fo,lam).real
    ns=index(SILICA,lam).real
    eps={**{k:(ne*ne,no*no,no*no) for k in LN_NAMES},'clad':((ns+sec.clad_dn)**2,)*3,'box':(ns*ns,)*3,
         'substrate':(index(SUBSTRATE_FILES[sec.substrate],lam).real**2,)*3,'air':(1.0,)*3}
    if sec.rail_t is not None:
        au=index(GOLD,lam)**2
        eps['rail_l']=eps['rail_r']=(au,)*3
    return eps


def cladding_index_max(sec,lam):
    """Largest refractive index among the non-LN dielectrics (TE0 must lie above it)."""
    return max(np.sqrt(v[0].real) for k,v in permittivity(sec,lam).items() if k not in LN_NAMES and not k.startswith('rail'))


def solve_modes(mesh,eps_by_region,lam,n_guess,num_modes=6,metallic=True,order=2,deps_xx=None,intorder=None):
    """femwell modes for a piecewise-constant diagonal permittivity. Returns (modes, seconds).

    deps_xx: optional per-element addition to eps_xx (length mesh.nelements; O02 convention B).
    intorder: quadrature order of the permittivity basis, which femwell's mode basis inherits. None keeps
    skfem's default for P0 (a 3-point rule), as used in O01; O02 passes 4 (as A17 / the femwell tutorials).
    """
    from femwell.maxwell.waveguide import compute_modes
    from skfem import Basis,ElementTriP0,ElementVector
    # femwell takes the diagonal-anisotropy branch when epsilon interpolates to (3, elements, points),
    # i.e. epsilon lives on a 3-component P0 basis; interior_dofs[k] holds component k of every element.
    basis0=Basis(mesh,ElementVector(ElementTriP0(),3),**({} if intorder is None else {'intorder':intorder}))
    lossy=any(np.iscomplexobj(np.asarray(v)) and np.any(np.imag(v)) for v in eps_by_region.values())
    eps=basis0.zeros(dtype=complex if lossy else float)
    for name,v in eps_by_region.items():
        if name not in mesh.subdomains:continue   # e.g. rail_l on a half domain
        elems=mesh.subdomains[name]
        for k in range(3):eps[basis0.interior_dofs[k,elems]]=v[k]
    if np.any(eps==0):raise ValueError('element without a material')
    if deps_xx is not None:
        eps=eps.astype(np.result_type(eps.dtype,np.asarray(deps_xx).dtype))
        eps[basis0.interior_dofs[0]]+=deps_xx
    t=time.perf_counter()
    modes=compute_modes(basis0,eps,wavelength=lam,num_modes=num_modes,order=order,metallic_boundaries=metallic,n_guess=n_guess)
    return modes,time.perf_counter()-t


def pick_te0(modes,n_floor,te_min=0.8):
    """O01 BOUNDARY 0.2: highest Re n_eff above n_floor with TE fraction >= te_min; None if absent."""
    ok=[m for m in modes if m.n_eff.real>n_floor and m.te_fraction>=te_min]
    return max(ok,key=lambda m:m.n_eff.real) if ok else None


def solve_te0(sec,mesh,lam,deps_xx=None,intorder=None):
    """TE0 of a Section at wavelength lam (um): n_eff, TE fraction, loss, and a summary of all modes found."""
    eps=permittivity(sec,lam)
    n_guess=float(np.sqrt(eps['ridge'][1]))   # shift k0^2 n_o^2 (boundary 0.2)
    modes,sec_s=solve_modes(mesh,eps,lam,n_guess,deps_xx=deps_xx,intorder=intorder)
    floor=cladding_index_max(sec,lam)
    te0=pick_te0(modes,floor)
    out={'lam_um':lam,'solve_s':sec_s,'n_floor':floor,'n_elements':int(mesh.nelements),
         'modes':[{'neff':complex(m.n_eff),'te_fraction':float(m.te_fraction)} for m in modes],'found':te0 is not None}
    if te0 is not None:
        n=complex(te0.n_eff)
        out.update(neff=n,te_fraction=float(te0.te_fraction),
                   loss_dB_per_cm=float(20/np.log(10)*2*np.pi/lam*abs(n.imag)*1e4),mode=te0)
    return out


def group_index(n_minus,n_0,n_plus,lam,delta):
    """ng = n - lam dn/dlam by central difference (lam, delta in the same unit)."""
    return n_0-lam*(n_plus-n_minus)/(2*delta)


# ---- G1 method gate: anisotropic slab with an exact dispersion relation --------------------------------

def slab_exact(eps_core,n_clad,d,lam,pol):
    """Fundamental mode of a symmetric slab (thickness d, core eps = (xx, yy, zz), isotropic cladding).

    TE0 (E along x) sees eps_xx only: kappa tan(kappa d/2) = gamma.
    TM0 (H along x): H'' term from d/dy[(1/eps_zz) dH/dy] + (k0^2 - beta^2/eps_yy) H = 0, so
    kappa = k0 sqrt(eps_zz/eps_yy) sqrt(eps_yy - n^2) and (kappa/eps_zz) tan(kappa d/2) = gamma/n_clad^2.
    """
    k0=2*np.pi/lam;exx,eyy,ezz=eps_core;nc2=n_clad**2
    if pol=='TE':
        top=exx;kap=lambda n:k0*np.sqrt(exx-n*n);w=lambda n:(kap(n),1.0,1.0)
    else:
        top=eyy;kap=lambda n:k0*np.sqrt(ezz/eyy)*np.sqrt(eyy-n*n);w=lambda n:(kap(n),ezz,nc2)
    def g(n):
        k,a,b=w(n);gam=k0*np.sqrt(n*n-nc2)
        return k/a*np.sin(k*d/2)-gam/b*np.cos(k*d/2)
    # Fundamental mode: kappa d/2 in (0, pi/2), where g > 0 at the low end (cos = 0 or gamma = 0) and
    # g = -gamma/b < 0 as kappa -> 0, with a single sign change in between.
    half_pi=top-(np.pi/(k0*d))**2*(1.0 if pol=='TE' else eyy/ezz)   # n^2 at kappa d/2 = pi/2
    lo=np.sqrt(max(nc2,half_pi))+1e-12;hi=np.sqrt(top)-1e-12
    if not g(lo)>0>g(hi):raise ValueError('no fundamental slab mode')
    return brentq(g,lo,hi,xtol=1e-15,rtol=1e-15)


def slab_geometry(d_core=0.6,clad=2.0,half_width=1.0):
    return OrderedDict(core=box(-half_width,-d_core/2,half_width,d_core/2),
                       cladding=box(-half_width,-d_core/2-clad,half_width,d_core/2+clad))


def build_slab_mesh(m,d_core=0.6,clad=2.0,half_width=1.0):
    from femwell.mesh import mesh_from_OrderedDict
    from skfem.io.meshio import from_meshio
    return from_meshio(mesh_from_OrderedDict(slab_geometry(d_core,clad,half_width),
                                             {'core':{'resolution':.02*m,'distance':.5}},default_resolution_max=.4*m))


def solve_slab(mesh,lam,pol,ln='CLN'):
    """FEM slab mode for the G1 gate: TE with PEC walls, TM with natural (PMC) walls."""
    fe,fo=LN_FILES[ln]
    ne,no=index(fe,lam).real,index(fo,lam).real;ns=index(SILICA,lam).real
    eps={'core':(ne*ne,no*no,no*no),'cladding':(ns*ns,)*3}
    modes,sec_s=solve_modes(mesh,eps,lam,n_guess=no,metallic=(pol=='TE'))
    want=(lambda m:m.te_fraction>0.5) if pol=='TE' else (lambda m:m.te_fraction<0.5)
    cand=[m for m in modes if want(m) and m.n_eff.real>ns]
    best=max(cand,key=lambda m:m.n_eff.real) if cand else None
    exact=slab_exact((ne*ne,no*no,no*no),ns,0.6,lam,pol)
    return {'neff':None if best is None else complex(best.n_eff),'te_fraction':None if best is None else float(best.te_fraction),
            'exact':exact,'solve_s':sec_s,'n_elements':int(mesh.nelements)}


