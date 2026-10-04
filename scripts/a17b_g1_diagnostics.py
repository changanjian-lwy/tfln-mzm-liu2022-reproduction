"""A17b G1 method diagnostics on the Tuncer CPW (see A17b BOUNDARY, preamble and section 8). Prints ratios only.

  shift  shift-invert sensitivity: for each order and frequency on the (4,1) mesh, the largest relative change of
         neff and |Z0| when the solve is repeated with the default shift or n_guess = (1, 1.02, 0.98) x n0 (n0 = the
         default solve); then the same at 50 MHz on (4,1) and (2,1). Low-frequency rows vary between runs: that
         drift is the finding.
  rate   mesh-halving changes at 5 GHz, d=1 (relative to the finer mesh): order 1 m=1..0.125, order 2 m=1..0.25.
  cost   one solve at 1 GHz: python scripts/a17b_g1_diagnostics.py cost M D ORDER (mesh size, DOFs, time, peak RSS).
"""
import resource,sys
from tfln_mzm.cpw_fem import CPW,build_mesh,solve

C=CPW(w_sig=7.0,gap=10.0,w_gnd=100.0,t_metal=0.8,sigma=6e7,eps_sub=13.0)
ch=lambda a,b:abs(a-b)/abs(b)


def shift_spread(mesh,f,order):
    d0=solve(C,mesh,f,order=order);n0=d0['neff'];out=[]
    for s in (None,1.0,1.02,0.98):
        r=solve(C,mesh,f,order=order,n_guess=None if s is None else s*n0)
        out.append((ch(r['neff'],n0),ch(r['Z0_pi'],d0['Z0_pi']),abs(r['gamma_rlgc']/r['gamma_fem']-1)))
    return [max(o[i] for o in out) for i in range(3)]


def shift():
    mesh=build_mesh(C,4.0,1)
    for order in (1,2):
        for f in (5e7,1.5e8,5e8,2e9,5e9):
            dn,dz,g=shift_spread(mesh,f,order)
            print(f'order {order} (4,1) {f/1e9:6.3f} GHz: max|dneff| {dn:.1e}  max|dZ0| {dz:.1e}  worst 2b gamma {g:.3f}',flush=True)
    for order in (1,2):
        for m in (4.0,2.0):
            dn,dz,g=shift_spread(build_mesh(C,m,1),5e7,order)
            print(f'order {order} ({m:g},1) 0.050 GHz: max|dneff| {dn:.1e}  max|dZ0| {dz:.1e}  worst 2b gamma {g:.3f}',flush=True)


def rate():
    for order,ms in ((1,(1,.5,.25,.125)),(2,(1,.5,.25))):
        prev=None
        for m in ms:
            r=solve(C,build_mesh(C,m,1),5e9,order=order)
            if prev:
                p=prev[1]
                print(f"order {order} m {prev[0]:g}->{m:g}: |Z0| {ch(abs(p['Z0_pi']),abs(r['Z0_pi'])):.3%}  nm {ch(p['nm'],r['nm']):.3%}  "
                      f"alpha {ch(p['alpha_Np_per_m'],r['alpha_Np_per_m']):.3%}  R {ch(p['R'],r['R']):.2%}  L {ch(p['L'],r['L']):.2%}  "
                      f"C {ch(p['C'],r['C']):.2%}",flush=True)
            prev=(m,r)


def cost(m,d,order):
    r=solve(C,build_mesh(C,m,d),1e9,order=order)
    print(f"(m,d)=({m:g},{d}) order {order}: {r['n_elements']} elements, {r['n_dofs']} DOFs, {r['solve_s']:.1f} s, "
          f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9:.2f} GB")


if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='cost':cost(float(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
    else:{'shift':shift,'rate':rate}[mode]()
