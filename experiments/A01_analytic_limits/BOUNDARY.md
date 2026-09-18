# A01 analytical limits (declared before execution)

Parent: A00 / commit 884a480. Scope: Eq.(1)–(4) and terminated-line S11 numerical implementation; not Figure 3 agreement.

Change: add modular analytical solver. No fitting and no changes to A00 missing-data declarations. Test scenarios are explicit numerical idealizations, not recovered author FEM data.

Fixed baseline: L=0.006 m, Vg=1 V, Zg=50 ohm, ng=2.25. Default ideal provider: Z0=50 ohm, nm=2.25, alpha=0 Np/m. The test matrix varies only the quantities needed by each declared limit: Z0={35,50,65}, ZL={0,20,40,50,80,open}, f=0 and 1 MHz–200 GHz, alpha={0,20,100} Np/m, mismatch=0.15, short-length limit L=1e-12 m. These are unit-test coordinates, not a paper parameter fit.

Success: B01–B09 in docs/BOUNDARIES.md pass at atol=1e-10 and rtol=1e-8 for non-singular points. Analytical complex Vavg must also match independent direct quadrature of printed Eq.(1), including loss and mismatch. S11 independently checked using ABCD input impedance.

Failure: any limit mismatch, nonpassive reflection, sign/unit error, silent zero-reference normalization, or inconsistent complex integral. Report the failing test; do not tune paper inputs.

Cannot prove: actual FEM parameters, Figure 3 numerical agreement, fabricated-device behavior or bandwidth extraction. B10–B12 are deferred to A02.
