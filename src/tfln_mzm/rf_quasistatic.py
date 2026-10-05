"""Quasi-static RF cross-section of the main CPW and the capacitively loaded line (SC-06, Track Q). Not imported by L0/L1.

Perfect conductors, 2D, high-frequency (external-inductance) limit (Q01/Q02 BOUNDARY). Half domain x >= 0 about the
signal centre: x = 0 is natural (even symmetry: signal at V, grounds at 0); the other outer edges are natural too
(zero net charge, the two-conductor line; Q01 BOUNDARY section 5 replaces the boundary's phi = 0 shield), and
the domain size is checked by doubling it. Every dof of a conductor region or conductor line is fixed.
C' = eps0 * phi^T A phi / V^2 for the half section, doubled for the full section; the 2D energy is scale-free, so
lengths in um give C' in F/m. L = 1 / (c^2 C0') with C0' the vacuum capacitance of the current-carrying conductors
(signal and grounds; T-heads and necks carry no longitudinal current, ref [10]). eps_LN = diag(eps33, eps11) in
(x, y) because x is crystal Z and y crystal X (optical_fem). Material constants: data/external/dc_materials and
rf_materials SOURCE.md.

Layer stack (y up, y = 0 at the top of the bonding SiO2, as in O01): quartz 500 um, bonding SiO2 2 um, LN (300 nm
etched slab, 600 nm ridge or unetched film), PECVD SiO2 100 nm, BCB of thickness t outside the window around the
arm, main electrodes on the BCB, T-heads on the PECVD top with a 3 um gap centred on the ridge.
"""
from collections import OrderedDict
from dataclasses import dataclass,replace
import numpy as np
from scipy.special import ellipk
from shapely.affinity import translate
from shapely.geometry import LineString,Polygon,box  # noqa: F401  (box is re-used by tests)
from shapely.ops import unary_union

EPS0=8.8541878128e-12
C_LIGHT=299792458.0
LN_EPS={'eps33':27.9,'eps11':44.0}            # clamped, dc_materials/SOURCE.md
SUBSTRATE_EPS={'quartz':4.5,'fused':3.8}
HEAD_GAP=3.0


# ---------------------------------------------------------------- exact conformal-mapping references (Q01 G1)
def k_cpw(s,g,wg=np.inf):
    """Modulus of a zero-thickness CPW with finite grounds (Ghione & Naldi, IEEE T-MTT 35(3), 1987)."""
    a,b=s/2,s/2+g
    if not np.isfinite(wg):return a/b
    c=b+wg
    return (a/b)*np.sqrt((1-b*b/(c*c))/(1-a*a/(c*c)))


def c_cpw_exact(s,g,wg=np.inf,eps_below=1.0,eps_above=1.0):
    """Full-section capacitance (F/m) of a zero-thickness CPW between two half-spaces: 2 eps0 (e1 + e2) K(k)/K(k').

    An anisotropic half-space diag(ex, ey) acts as an isotropic one with sqrt(ex ey) (scale y by sqrt(ex/ey))."""
    k=k_cpw(s,g,wg)
    return 2*EPS0*(eps_below+eps_above)*ellipk(k*k)/ellipk(1-k*k)


# ---------------------------------------------------------------- meshing
def _vertices(geom):
    """(vertex, unit vectors to its neighbours along the outline) for a polygon's corners or a line's two ends."""
    if isinstance(geom,Polygon):
        p=np.asarray(geom.exterior.coords)[:-1];n=len(p)
        return [(p[i],[_unit(p[i-1]-p[i]),_unit(p[(i+1)%n]-p[i])]) for i in range(n)]
    for g in getattr(geom,'geoms',[geom]):
        p=np.asarray(g.coords)
        return [(p[0],[_unit(p[1]-p[0])]),(p[-1],[_unit(p[-2]-p[-1])])]


def _unit(v):return v/np.linalg.norm(v)


def build_mesh(shapes,conductors,m,X,thin=(),corner_res=0.02,face_res=0.5,far=None):
    """Mesh an OrderedDict of shapes (earlier entries win) with graded refinement at conductor corners and faces.

    Each conductor outline gets face lines (size face_res m, growing to 5 m within 30 um) and every corner off the
    x = 0 symmetry plane a short mark along each adjacent edge (size corner_res m, growing to face_res m within
    4 um). thin: (name, thickness) layers whose size is capped at min(20 thickness, 5) m. The far field is capped
    at `far` m (default X/40). Returns (mesh, {conductor: [mark and face line names lying on it]})."""
    from femwell.mesh import mesh_from_OrderedDict
    from skfem.io.meshio import from_meshio
    far=X/40 if far is None else far
    ell=corner_res*m/2
    res={};lines=OrderedDict();owned={}
    for name in conductors:
        geom=shapes[name];owned[name]=[]
        for i,(p,dirs) in enumerate(_vertices(geom)):
            if abs(p[0])<1e-9:continue                             # symmetry plane: not an edge
            for j,u in enumerate(dirs):
                key=f'{name}__c{i}_{j}';lines[key]=LineString([tuple(p),tuple(p+ell*u)])
                res[key]={'resolution':corner_res*m,'distance':4.0,'SizeMax':face_res*m};owned[name].append(key)
        if isinstance(geom,Polygon):
            key=f'{name}__face';lines[key]=LineString(list(geom.exterior.coords))
            res[key]={'resolution':face_res*m,'distance':30.0,'SizeMax':5*m}
            res[name]={'resolution':5*m,'distance':1e-3}
        else:
            res[name]={'resolution':face_res*m,'distance':30.0,'SizeMax':5*m}
    for name,th in thin:
        res[name]={'resolution':min(20*th,5.0)*m,'distance':1e-3}
    allshapes=OrderedDict(list(lines.items())+list(shapes.items()))
    raw=mesh_from_OrderedDict(allshapes,res,default_resolution_max=far*m)
    return from_meshio(raw),owned


# ---------------------------------------------------------------- solve
def solve_capacitance(mesh,eps,fixed_regions,fixed_lines=(),extra_zero_facets=None,outer='natural'):
    """Half-section capacitance (F/m) of the conductor set at potential 1 against everything at 0.

    eps: {region: (exx, eyy)} for every dielectric region; conductor regions get 1 (their dofs are all fixed).
    fixed_regions / fixed_lines: [(name, potential)]. outer='natural' (zero normal D: zero net charge, the
    two-conductor line) or 'zero' (phi = 0 on the outer edges except x = 0, a grounded shield)."""
    from skfem import Basis,BilinearForm,ElementTriP0,ElementTriP2,condense,solve
    basis=Basis(mesh,ElementTriP2());b0=basis.with_element(ElementTriP0())
    exx=np.full(b0.N,np.nan);eyy=np.full(b0.N,np.nan)
    for name,_ in fixed_regions:
        d=b0.get_dofs(elements=name).all();exx[d]=1.0;eyy[d]=1.0
    for name,(ex,ey) in eps.items():
        if name not in mesh.subdomains:continue
        d=b0.get_dofs(elements=name).all();exx[d]=ex;eyy[d]=ey
    if np.isnan(exx).any():raise ValueError(f'{int(np.isnan(exx).sum())} elements without a permittivity')

    @BilinearForm
    def a(u,v,w):return w.exx*u.grad[0]*v.grad[0]+w.eyy*u.grad[1]*v.grad[1]
    A=a.assemble(basis,exx=b0.interpolate(exx),eyy=b0.interpolate(eyy))
    x=basis.zeros();D=[]
    sets=[]
    if outer=='zero':
        sets.append((basis.get_dofs(facets=mesh.facets_satisfying(lambda p:p[0]>1e-9,boundaries_only=True)).all(),0.0))
    elif outer!='natural':raise ValueError(outer)
    for name,v in fixed_regions:sets.append((basis.get_dofs(elements=name).all(),v))
    for name,v in fixed_lines:
        if name in mesh.boundaries:sets.append((basis.get_dofs(facets=mesh.boundaries[name]).all(),v))
    if extra_zero_facets is not None:sets.append((basis.get_dofs(facets=extra_zero_facets).all(),0.0))
    for dofs,v in sets:
        dofs=np.asarray(dofs);x[dofs]=v;D.append(dofs)
    D=np.unique(np.concatenate(D))
    phi=solve(*condense(A,x=x,D=D))
    return EPS0*float(phi@(A@phi)),basis,phi


# ---------------------------------------------------------------- Q01: main CPW in vacuum
def main_cpw_shapes(s=80.0,g=20.0,wg=150.0,t=4.0,X=2000.0,y0=0.0,thin=False):
    """Signal and one ground (half section) in vacuum; thin=True gives zero-thickness lines at y0."""
    shapes=OrderedDict()
    if thin:
        shapes['sig']=LineString([(0,y0),(s/2,y0)]);shapes['gnd']=LineString([(s/2+g,y0),(s/2+g+wg,y0)])
    else:
        shapes['sig']=box(0,y0,s/2,y0+t);shapes['gnd']=box(s/2+g,y0,s/2+g+wg,y0+t)
    shapes['air']=box(0,-X,X,X)
    return shapes


def half_space_shapes(s=80.0,g=20.0,wg=150.0,X=2000.0):
    """Zero-thickness CPW on the interface y = 0 between 'above' (y > 0) and 'below' (y < 0)."""
    shapes=OrderedDict()
    shapes['sig']=LineString([(0,0),(s/2,0)]);shapes['gnd']=LineString([(s/2+g,0),(s/2+g+wg,0)])
    shapes['above']=box(0,0,X,X);shapes['below']=box(0,-X,X,0)
    return shapes


def inductance(c0_full):
    """External inductance (H/m) from the full-section vacuum capacitance."""
    return 1.0/(C_LIGHT**2*c0_full)


# ---------------------------------------------------------------- Q02: Liu cross-section
@dataclass(frozen=True)
class RFSection:
    """Q02 BOUNDARY section 1 table (lengths um). kind: 'U' unloaded, 'H' through T-heads, 'N' through necks."""
    t_bcb:float=1.5
    kind:str='H'
    w_t:float=2.0
    h_t:float=0.5
    window:float|None=None        # BCB window half-width about the ridge centre; None = aligned with electrodes
    etch:str='full'               # 'full' (300 nm slab everywhere) or 'gap' (600 nm outside the electrode gap)
    wg:float=150.0
    s:float=80.0
    g:float=20.0
    t_au:float=4.0
    t_film:float=0.6
    etch_depth:float=0.3
    t_clad:float=0.1
    t_box:float=2.0
    t_sub:float=500.0
    ridge_top:float=1.0
    theta_deg:float=75.0

    @property
    def xa(s):return s.s/2+s.g/2
    @property
    def W(s):return s.g/2 if s.window is None else s.window


def liu_shapes(sec,X):
    """Ordered shapes (conductors first) and the list of (thin layer, thickness) for build_mesh."""
    xa=sec.xa;ys=sec.t_film-sec.etch_depth
    run=sec.etch_depth/np.tan(np.radians(sec.theta_deg))
    ridge=Polygon([(xa-sec.ridge_top/2-run,ys),(xa+sec.ridge_top/2+run,ys),(xa+sec.ridge_top/2,sec.t_film),
                   (xa-sec.ridge_top/2,sec.t_film)])
    x1,x2=sec.s/2,sec.s/2+sec.g
    if sec.etch=='full':ln=unary_union([box(0,0,X,ys),ridge])
    elif sec.etch=='gap':ln=unary_union([box(0,0,x1,sec.t_film),box(x1,0,x2,ys),box(x2,0,X,sec.t_film),ridge])
    else:raise ValueError(sec.etch)
    top=box(0,ys,X,X)
    clad=ln.buffer(sec.t_clad,join_style='mitre').intersection(top).intersection(box(0,0,X,X)).difference(ln)
    below=unary_union([ln,clad,box(0,-sec.t_box-sec.t_sub,X,0)])
    y_gap=ys+sec.t_clad                                   # PECVD top inside the electrode gap
    y_el=(ys if sec.etch=='full' else sec.t_film)+sec.t_clad+sec.t_bcb   # electrode bottom
    win=box(xa-sec.W,-1,xa+sec.W,X)
    shapes=OrderedDict()
    shapes['sig']=box(0,y_el,x1,y_el+sec.t_au);shapes['gnd']=box(x2,y_el,x2+sec.wg,y_el+sec.t_au)
    cond=['sig','gnd']
    if sec.kind in ('H','N'):
        hs=box(xa-HEAD_GAP/2-sec.w_t,y_gap,xa-HEAD_GAP/2,y_gap+sec.h_t)
        hg=box(xa+HEAD_GAP/2,y_gap,xa+HEAD_GAP/2+sec.w_t,y_gap+sec.h_t)
        if sec.kind=='N':
            hs=unary_union([hs,*_neck(sec,y_gap,y_el,-1)]);hg=unary_union([hg,*_neck(sec,y_gap,y_el,+1)])
        shapes['head_s']=hs;shapes['head_g']=hg;cond+=['head_s','head_g']
    shapes['ridge']=ridge.difference(box(-1,-1,0,X))
    shapes['ln']=ln
    shapes['clad']=clad
    if sec.t_bcb>0:
        bcb=translate(below,yoff=sec.t_bcb).intersection(box(0,0,X,X)).difference(below).difference(win)
        shapes['bcb']=bcb
    shapes['box']=box(0,-sec.t_box,X,0)
    shapes['quartz']=box(0,-sec.t_box-sec.t_sub,X,-sec.t_box)
    shapes['air']=box(0,-X,X,X)
    thin=[('ln',ys),('clad',sec.t_clad),('ridge',ys)]
    if sec.t_bcb>0:thin.append(('bcb',sec.t_bcb))
    return shapes,cond,thin


def _neck(sec,y_gap,y_el,side):
    """Neck conductor pieces from the head to the main electrode (side -1: signal, +1: ground)."""
    xa=sec.xa;h=sec.h_t
    xe=sec.s/2 if side<0 else sec.s/2+sec.g              # electrode inner edge
    xw=xa+side*sec.W                                     # BCB wall
    xh=xa+side*(HEAD_GAP/2+sec.w_t)                      # head outer edge
    lo,hi=sorted((xw,xh))
    pieces=[box(lo,y_gap,hi,y_gap+h)]
    if sec.t_bcb>0 or abs(y_el-y_gap)>1e-9:
        ledge_top=y_gap+sec.t_bcb
        wl,wh=sorted((xw,xw-side*h))                     # strip on the gap side of the wall
        pieces.append(box(wl,y_gap,wh,max(ledge_top,y_gap)+h))
        if abs(xw-xe)>1e-9:
            a,b=sorted((xe,xw-side*h))
            pieces.append(box(a,ledge_top,b,ledge_top+h))
            c,d=sorted((xe,xe-side*h))
            pieces.append(box(c,ledge_top,d,y_el))      # riser to the electrode bottom
    return pieces


def liu_eps(eps_bcb=2.65,ln_scale=1.0,eps_sio2=3.9,substrate='quartz',vacuum=False):
    """{region: (exx, eyy)}; vacuum=True sets every dielectric to 1 (for L)."""
    if vacuum:
        return {k:(1.0,1.0) for k in ('ridge','ln','clad','bcb','box','quartz','air')}
    ln=(LN_EPS['eps33']*ln_scale,LN_EPS['eps11']*ln_scale)
    q=SUBSTRATE_EPS[substrate]
    return {'ridge':ln,'ln':ln,'clad':(eps_sio2,)*2,'bcb':(eps_bcb,)*2,'box':(eps_sio2,)*2,'quartz':(q,q),
            'air':(1.0,1.0)}


def loaded_c(cu,ch,cn,duty=0.9,fn=0.05):
    """C_avg = (D - f_n) C_H + f_n C_N + (1 - D) C_U (Q02 BOUNDARY)."""
    return (duty-fn)*ch+fn*cn+(1-duty)*cu
