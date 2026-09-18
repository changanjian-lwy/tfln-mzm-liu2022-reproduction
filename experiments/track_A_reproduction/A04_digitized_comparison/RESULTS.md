# A04 results

Status: PARTIAL_COMPARISON_FAIL. Full Figure 3 reproduction is not established. No input parameter was fitted to target curves.

| Panel | Load | Scored points | MAE | Within readout allowance | Diagnostic pass |
|---|---|---:|---:|---:|---|
| fig3a_inset | 20 | 392 | 0.0088 V | 74.7% | False |
| fig3a_inset | 40 | 360 | 0.0041 V | 95.0% | True |
| fig3a_inset | 50 | 130 | 0.0065 V | 79.2% | False |
| fig3a_inset | 80 | 310 | 0.0128 V | 46.5% | False |
| fig3a_inset | open | 364 | 0.0339 V | 17.0% | False |
| fig3b | 20 | 733 | 0.7572 dB | 1.8% | False |
| fig3b | 40 | 787 | 0.0677 dB | 98.0% | True |
| fig3b | 50 | 784 | 0.0703 dB | 98.7% | True |
| fig3b | 80 | 787 | 0.1223 dB | 85.3% | False |
| fig3b | open | 787 | 0.4274 dB | 30.2% | False |
| fig3c | 20 | 716 | 2.2955 dB | 24.7% | False |
| fig3c | 40 | 760 | 8.4283 dB | 0.0% | False |
| fig3c | 50 | 718 | 5.3841 dB | 22.4% | False |
| fig3c | 80 | 760 | 6.3927 dB | 7.2% | False |

Figure 3(b) 40-ohm and 50-ohm responses meet the predeclared diagnostic on retained points; Figure 3(a) inset 40-ohm also meets it. The other series fail. This is not a full-band or complete-panel validation: large source gaps and target legend/axis contamination are excluded.

S11 differences are substantial for all four loads. The current representation uses real-valued Z0 read from an unlabeled scalar impedance plot and constant nm=2.25; it does not recover complex impedance, nm(f), feedline behavior or the authors’ precise reference-plane and low-frequency normalization choices. These are candidate explanations, not proven causes. Passing ABCD and passivity tests does not prove the chosen input provider matches the authors’ device.

The normalized 20-ohm curve shows an approximately systematic offset; do not fix it by scaling the curve. Next audit must distinguish source normalization, reference frequency, impedance interpretation and missing dispersion. The low-frequency provider is explicitly assumed, not measured or extracted.

A03 axis-tick QA correction occurred after the initial comparison; results_before_axis_tick_QA.json retains that initial metric report. No acceptance or physical parameters changed. Updated comparison image was inspected.

All 27 implementation tests passed. Source-data CSV hash, all residual rows, exclusion counts and per-series diagnostics are archived. No bandwidth is claimed for the incomplete digitized provider. Next step: a separately bounded convention/input audit, followed by full Figure 3(a) extraction if needed.
