# Track O: optical side from the cross-section (SC-05)

Authorized 2026-10-05 12:07 (decision D15). Every model input used elsewhere in this project (α, nm, Z0) comes from the paper's own simulated curves; Track O instead computes optical quantities from the cross-section the paper describes: O01 the TE0 group index (paper: ng ∼ 2.25), O02 the half-wave-voltage–length product (paper: Vπ = 2.7 V over 6 mm, VπL = 1.62 V·cm, open terminal, 1 MHz triangular wave). Nothing here enters a Figure 3 comparison, and no RF FEM is involved (SC-01 stays paused).

What the paper leaves out (wavelength, sidewall angle, which width is 1 µm, T-rail thickness, duty cycle and position) is handled as declared ranges, never fitted. Where the cited predecessor paper [10] (Liu et al., Chin. Opt. Lett. 19, 060016, 2021) states a value, it is used as `external_reference` and labelled as such. Material data are in [`data/external/optical_materials/`](../../data/external/optical_materials/SOURCE.md). Results are intervals, and verdicts are only "compatible", "incompatible" or "undetermined".

Each experiment follows [EXPERIMENT_TEMPLATE.md](../EXPERIMENT_TEMPLATE.md); BOUNDARY.md is committed before the script, and the script before the run. FEM runs obey `BUDGET.json` (BOUNDARIES §10).

| Experiment | Grade | Claim boundary |
|---|---|---|
| [O01 optical TE0 group index](O01_optical_mode_ng/RESULTS.md) | DIAGNOSTIC_ONLY | Undetermined by the pre-declared rule. TE0 group index 2.2814–2.2861 over 13 configurations (wavelength 1530–1570 nm, sidewall 60–90°, 1 µm as top or bottom width, PECVD index +0–0.02, two substrates, T-rails absent or 0.2–0.8 µm thick); numerical uncertainty 7.8e-6. The paper's "~2.25" lies 0.026–0.036 below, more than ten times any declared factor; 5% MgO branch (not in the interval) 2.2739 |
