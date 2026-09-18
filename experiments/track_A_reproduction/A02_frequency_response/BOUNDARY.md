# A02: frequency response and bandwidth extraction

Parent: A01 analytical limits. Claim: QUALITATIVE_ONLY, not Figure 3 numerical reproduction.

Only swept physical variable: ZL=20/40/50/80 ohm/open. Fixed paper values: Vg=1 V, Zg=50 ohm, L=6 mm. Explicit assumed provider: constant Z0=50 ohm, nm=ng=2.25, alpha_dB_per_mm=0.05*sqrt(f_GHz). The coefficient is an illustrative choice, not a digitized or fitted paper value. No feedline, impedance dispersion, dielectric-loss term or measured inputs.

Reference f0=1 MHz for each load. Range f0 to 200 GHz, coarse step <=0.1 GHz, fine <=0.05 GHz. Bandwidth is the first downward -3 dB crossing relative to each load's own f0 response, not relative to peak. Interpolate finite dB brackets; multiple crossings listed. A sampled exact threshold counts as first threshold contact; a tangency is labelled as contact, not asserted to cross. If a bracket includes -infinity return a bracket-only result requiring refinement; do not invent an interpolated frequency. If no threshold contact/crossing exists, return null bandwidth, censored=true and lower bound fmax (conditional on grid resolution).

Success: B10 grid refinement differences <=0.1 GHz for resolved crossings; B11 no-crossing censor; B12 first crossing despite later recrossings. Independent mismatch-only sinc first-lobe root validates bandwidth extraction. Coarse/fine censor status and event counts must agree. A02 figures preserve absolute Vavg and normalized response separately. Zero S11 stays -infinity in data; a plotting note explains the invisible matched curve.

Failure: invented crossing, changing normalization to improve bandwidth, hiding multiple crossings, or convergence failure. Cannot prove author FEM inputs, Figure 3 pointwise match, measured flatness or hardware bandwidth.
