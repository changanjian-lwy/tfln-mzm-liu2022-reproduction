# A02 results

Grade: QUALITATIVE_ONLY. Simulation completed with the predeclared assumed provider. No fitting or parameter adjustment was performed. Figure 3 quantitative validation remains pending.

All 25 tests pass, including five new bandwidth tests covering B10–B12, invalid grids, threshold contacts and infinite null brackets. The mismatch-only first-lobe sinc root independently agrees with the extractor. Coarse/fine grids have consistent event counts and censor status for all five loads.

| Load | First -3 dB bandwidth, assumed model | Grid refinement shift |
|---|---|---|
| 20 ohm | >200 GHz, censored | No crossing on either grid |
| 40 ohm | >200 GHz, censored | No crossing on either grid |
| 50 ohm | >200 GHz, censored | No crossing on either grid |
| 80 ohm | 63.4627 GHz | 0.000846 GHz |
| open | 4.1121 GHz | 0.0000541 GHz |

These values describe the chosen illustrative loss and constant-impedance model, not the paper's device. The lower-load normalized response rises initially while its low-frequency absolute voltage is smaller. Constant Z0=Zg causes matched-load S11 to be exactly zero at all frequencies and eliminates the impedance-dispersion features in paper Fig.3(c). Thus this is an explicit example of the limitation of the assumed provider, not successful reproduction of that panel. Open-load S11 is additionally plotted as a model diagnostic, not as a Figure 3(c) target.

`curves.csv` preserves complex Vavg, magnitude, normalized EO response and S11 in dB for 4001 frequencies per load. Matched-load S11 is stored as -inf rather than an arbitrary floor. `results.json` includes the actual parameters, reference frequency, all detected events, brackets, censor flags and refinement errors. `response.png` is labelled as assumed-input simulation and was visually inspected.

Reproduce with `python scripts/run_a02.py`; run checks with `python -m unittest discover -s tests -v`. Censored bounds remain conditional on frequency resolution; this convergence check is not a proof that arbitrarily narrow unsampled features cannot exist.

Next: A03 digital extraction of Figure 2 input curves and Figure 3 reference curves, with source calibration and uncertainty recorded before model comparison. A00's missing original FEM arrays remain unresolved. Do not tune this assumed provider to the reference points and call it independent validation.
