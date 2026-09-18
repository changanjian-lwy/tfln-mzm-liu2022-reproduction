# A05: narrowed causes, not a solved reproduction

Grade: DIAGNOSTIC_ONLY. A04 and its fail verdict are unchanged. No physical baseline was fitted or replaced.

## 1. Cross-panel reference is not exactly the assumed DC divider

Using author Fig.3(a) inset and Fig.3(b), independently of our voltage model, infer Vref=Vabsolute/10^(M/20). Values below are medians of raster-derived overlapping samples, not recovered author settings or experimental voltages.

| Load | Ideal divider V | Implied reference V | Equivalent offset dB | Original EO MAE dB | Output-calibrated diagnostic MAE dB |
|---|---:|---:|---:|---:|---:|
| 20 | 0.28571 | 0.30585 | +0.591 | 0.757 | 0.216 |
| 40 | 0.44444 | 0.44494 | +0.010 | 0.068 | 0.063 |
| 50 | 0.50000 | 0.49287 | -0.125 | 0.070 | 0.084 |
| 80 | 0.61538 | 0.60088 | -0.207 | 0.122 | 0.171 |
| open | 1.00000 | 0.94993 | -0.446 | 0.427 | 0.325 |

The 20-ohm systematic error is substantially reduced in this diagnostic, so reference normalization is a meaningful contributor. It does not explain all residuals or improve every load. Do not promote the inferred voltages to paper-explicit inputs. Adjacent raster samples are correlated; the sample count is not the number of independent measurements. Approximate vertical uncertainty and 5–95% spread are in JSON; horizontal readout uncertainty is not included in this reference calculation.

The visible curves begin at about 1 GHz. A finite-frequency/first-sample reference is a candidate explanation, but neither the exact reference frequency nor the authors' implementation is identified here. The paper's general statement that f0 is usually DC or MHz must not be overwritten with an invented 1-GHz setting.

## 2. S11 cannot generally be fixed by phase alone

For each supplied real Z0 and attenuation, analytically bound S11 magnitude over every possible propagation phase. Then separately widen the bound using a sampled local source uncertainty envelope.

| Load | Supported points | Outside nominal bounds | Outside sampled source-uncertainty envelope |
|---|---:|---:|---:|
| 20 | 716 | 179 | 37 |
| 40 | 760 | 752 | 652 |
| 50 | 724 | 275 | 58 |
| 80 | 764 | 437 | 329 |

Conclusion: with the current real-Z0, resistive-load, uniform-line input representation, varying only microwave phase cannot reproduce many S11 ordinates. The sampled uncertainty envelope is a sensitivity screen, not a certified interval or a complete uncertainty model. It nevertheless makes a simple small readout-error explanation less plausible for the 40/80-ohm discrepancies.

## 3. Refine the reference-plane explanation

A matched lossless line extension multiplies S11 by a unit phasor and cannot change |S11|. Thus a mere phase reference-plane shift cannot fix magnitude errors. A lossy or mismatched feedline/network can change magnitude, but needs source data and an explicit new module; it cannot be inserted as an unexplained correction.

## Remaining unknowns and next step

Audit whether the plotted Z0 is a magnitude, real part or an effective quantity and whether the loss curve and Fig.3 S11 share a de-embedded reference plane. Full complex Z0, microwave dispersion, feedline S parameters and exact per-load normalization settings are not recoverable from the current readout. These missing inputs are hypotheses, not proof of author mistakes. A04 also has incomplete frequency coverage and only Fig.3(a)'s inset.

Next work should be a source/convention clarification or separately labelled sensitivity study for the missing complex/network inputs. Do not change the load labels or fit source impedance just to match the curves. Full reproduction remains unresolved.

Verification: 29 tests pass. New audit tests compare the analytic bounds against an independent 20001-phase sweep over 36 impedance/load/loss combinations, and verify lossless reference-plane magnitude invariance. `python scripts/run_a05.py` regenerates diagnostics. All original A04 outputs remain intact.
