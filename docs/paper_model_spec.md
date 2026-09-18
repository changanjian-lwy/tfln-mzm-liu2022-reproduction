# Traveling-wave reproduction specification (main track, RF implementation pending)

Primary source: Liu et al., DOI [10.1109/LPT.2022.3178214](https://doi.org/10.1109/LPT.2022.3178214), PDF p.2 (journal p.855), equations (1)–(4); Figure 3 on PDF p.3 (p.856). Equations below were visually checked against the local PDF. Optical propagation subscript is written consistently as `o` here.

\[
V_{avg}(\omega)=\frac1L\int_0^L
\frac{V_g(1+\rho_1)e^{i\beta_oL}
\left[e^{i(\beta_e-\beta_o)x}+\rho_2e^{-i(\beta_e+\beta_o)x}\right]}
{2\left[e^{i\beta_eL}+\rho_1\rho_2e^{-i\beta_eL}\right]}dx.\tag{1}
\]
\[
\rho_1=\frac{Z_0-Z_g}{Z_0+Z_g},\tag{2}\qquad
\rho_2=\frac{Z_L-Z_0}{Z_L+Z_0}.\tag{3}
\]

For implementation, plot the magnitude response:

\[
M(f)=20\log_{10}\left|\frac{V_{avg}(2\pi f)}{V_{avg}(2\pi f_0)}\right|.\tag{4, magnitude interpretation}
\]

The printed Eq.(4) omits explicit modulus bars. A complex phasor ratio must be converted to magnitude for a real dB curve; log of the complex number is not the displayed observable. Retain phase separately if needed.

## Conventions to verify before coding

- Vg is the Thevenin source amplitude before the 50 Ω source impedance. Day 1 uses the corresponding ideal DC divider, not an RF calculation.
- βo=ω n_g/c describes modulation-envelope walk-off, not the optical carrier wavenumber 2π n_e/λ0.
- The coordinate convention of Eq.(1) must be reconciled with a port-based forward-wave model. For the printed exponential form, choosing βe=ω n_m/c−iα yields attenuation for e^(−iβe L), α≥0. Check the passive matched limit rather than guessing the loss sign.
- ρ1 in the paper is the input-interface coefficient with the paper's sign convention. A wave incident on the generator from inside the line sees Γg=(Zg−Z0)/(Zg+Z0)=−ρ1. Do not substitute one for the other silently.
- Use a stable analytic exponential integral around zero argument (expm1 or series), then cross-check against independent numerical quadrature, including perfect velocity match.
- Eq.(1) at zero frequency and zero line loss must reduce to Vg ZL/(Zg+ZL), even for Z0≠Zg. The short-circuit DC normalization is zero, so a short load cannot be normalized to DC without a different definition.
- S11 is a supplementary transmission-line derivation; it is **not** a fifth numbered paper equation. With 50 Ω reference and γ=α+iω n_m/c, derive input reflection from the terminated line and verify both impedance-transform and multiple-reflection expressions agree. Keep reference plane, feedline omission, and open-load limit explicit.
- Use f0=1 MHz as an explicit project choice unless the comparison needs a different documented baseline. Every load uses its own Vavg(f0); also retain absolute volts so reduced low-frequency drive is visible.

## Figure 3 validation plan

| Panel | Targets | Validation |
|---|---|---|
| (a) average voltage magnitude | 20/40/50/80 Ω and open; 0–200 GHz, inset 0–50 GHz | DC ordering, initial behavior, convergence trend, extracted anchor points |
| (b) normalized EO response | same five loads | lower loads can give an initial relative rise; high-load/open response falls; ripple and high-frequency slope |
| (c) S11 | 20/40/50/80 Ω | magnitude in dB at source reference plane; compare level and ripple with uncertainty |

Stage A: constant Z0, matched velocity and explicitly assumed loss law demonstrate mechanism. This is qualitative validation only.

Stage B: digitize relevant Figure 2 inputs and Figure 3 outputs separately. Freeze the input parameters before comparing held-out output samples. Report mean absolute error, maximum error and physically meaningful feature errors against extraction uncertainty. Treat deep S11 nulls cautiously because a tiny frequency shift makes dB error large. Do not silently tune unknown inputs to every output point.

No precise Figure 3 numerical validation is claimed today. The 20/40/50/80 Ω legends are visually verified, but feature values and crossing frequencies must be digitized before entering an acceptance table. Author FEM dispersion and feedline details limit exact reproduction.

## Limited extension

Proposed study grid (all study choices, not author data): L=3, 6, 9, 12 mm; RF attenuation multiplier=0.5, 1, 1.5, 2; Δn=n_m−n_g=−0.1, −0.05, 0, 0.05, 0.1. Start with one-factor-at-a-time plots at fixed ZL, then a small combined grid only if useful. Keeping Z0 and overlap constant while varying length/loss/mismatch isolates mechanisms but does not represent a self-consistent geometry redesign.

Define bandwidth as the **first downward crossing of −3 dB relative to f0**, interpolated between bracketing samples and verified with a finer grid. Peaking and ripple can create multiple crossings: report this convention and optionally the later crossings; never choose the last crossing solely for a larger number. If there is no crossing by fmax, report `> fmax` with a censored flag, not `fmax` as a measured bandwidth. No extrapolation beyond digitized input coverage without a separate assumption label.
