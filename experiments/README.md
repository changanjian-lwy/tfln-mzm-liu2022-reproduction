# Experiment index

Read [PROJECT_CHARTER.md](../PROJECT_CHARTER.md) and [docs/BOUNDARIES.md](../docs/BOUNDARIES.md) (v2.0) before running a new experiment. From A06 on, every experiment uses [EXPERIMENT_TEMPLATE.md](EXPERIMENT_TEMPLATE.md), commits its BOUNDARY.md before running, and passes gates G0–G4. Every experiment owns BOUNDARY.md, INPUTS.md and RESULTS.md; preserve both successful and failed conclusions. A new physical assumption, provider or conflict branch gets a new experiment rather than overwriting a completed one.

Evidence check (gate G1): `python scripts/check_evidence.py --replay` regenerates every replayable output in a temporary copy and requires byte-identical files. Add `--images ../tmp/pdfs/a03` to include the A03 digitization, which needs local paper rasters that are not in Git.

| Track | Experiment | Layer | Grade | P1 review (2026-10-04) | Claim boundary |
|---|---|---|---|---|---|
| A | [A00 baseline](track_A_reproduction/A00_paper_baseline/RESULTS.md) | — | BLOCKED_BY_MISSING_DATA | Citable; replay identical | Preflight only |
| A | [A01 analytical limits](track_A_reproduction/A01_analytic_limits/RESULTS.md) | L0 | ANALYTIC_LIMIT_PASS | Citable; replay identical | Eq.(1)–(4) / S11 analytical checks, not Figure 3 agreement |
| A | [A02 frequency response](track_A_reproduction/A02_frequency_response/RESULTS.md) | L0 | QUALITATIVE_ONLY; B10–B12 pass | Citable for bandwidth logic; replay identical | Assumed provider, first crossing, censoring and convergence |
| A | [A03 digitization](track_A_reproduction/A03_digitization/RESULTS.md) | data | PARTIAL_EXTRACTION_COMPLETE | Citable; replay identical with local rasters | Raster curves with explicit gaps; Fig.2(b) and full Fig.3(a) missing |
| A | [A04 independent comparison](track_A_reproduction/A04_digitized_comparison/RESULTS.md) | L1 | PARTIAL_COMPARISON_FAIL | Citable; replay identical | Some EO curves agree; S11 unresolved |
| A | [A05 convention audit](track_A_reproduction/A05_convention_audit/RESULTS.md) | L1 | DIAGNOSTIC_ONLY | Citable; replay identical | Normalization contributor identified; phase-only S11 correction ruled out for many samples |
| A | [A06 Fig.2(b) + Fig.3(a) main digitization](track_A_reproduction/A06_fig2b_fig3a_main_digitization/RESULTS.md) | data | PARTIAL_EXTRACTION_COMPLETE | First experiment under v2 rules; boundary committed first (`b84ed7a`) | nm(100 GHz)=2.2500 and 50-ohm bandwidth 202.6 GHz at BCB 1.5 um; new conflict C-06 |
| A | [A07 C-06 bandwidth definition](track_A_reproduction/A07_c06_bandwidth_definition/RESULTS.md) | L1 | DIAGNOSTIC_ONLY | Boundary committed first (`42f09e0`) | 50 ohm first −3 dB crossing: model 141.82 vs Fig.3(b) 141.78 GHz; ripple-insensitive definitions are compatible with Fig.2(b)'s 202.6 GHz |
| A | [A08 C-06 matched provider](track_A_reproduction/A08_c06_matched_provider/RESULTS.md) | L1 + C-06b | DIAGNOSTIC_ONLY | Boundary committed first (`f714833`) | A matched 50-ohm line with Fig.2(a) loss gives −2.3 to −2.5 dB at the Fig.2(b) bandwidths, so it cannot explain that curve |
| A | [A09 qualitative statements](track_A_reproduction/A09_qualitative_statements/RESULTS.md) | L1 | QUALITATIVE_ONLY | Boundary committed first (`c4445e9`) | Initial-rise/drop statement reproduced; roll-off similarity undecidable for the model (Z0 gap at 10 GHz); the S11 statement fails on the author's own figure as written; new conflict C-07 (possible 40/80-ohm legend swap in Fig.3(c)) |
| B | Termination, length, RF loss and mismatch scans | L0/L1 | Planned | — | Separate from paper baseline |

P1 finding that applies to A00–A05: each experiment's BOUNDARY.md and RESULTS.md entered Git in the same commit, so "criteria fixed before results" rests on the documents' own statements, not on Git history. These experiments predate the v2 rules and are not rewritten. Their INPUTS.md files are retrospective records written on 2026-10-04.

## P2 queue (within the authorized scope, no FEM)

Order follows the minimum evidence list in [SOURCE_COVERAGE_MATRIX.md](../docs/SOURCE_COVERAGE_MATRIX.md) §4. Numbers are assigned when each BOUNDARY.md is committed.

1. ~~Digitize Fig.2(b) and the full Fig.3(a)~~ — done in A06.
1b. ~~C-06 definition diagnostic~~ — A07. 1c. ~~C-06 matched-provider diagnostic~~ — A08 (ruled out under Fig.2(a) loss).
2. ~~Qualitative statements~~ — A09. 2b. C-07 legend-swap diagnostic for Fig.3(c) (A10).
3. C-02 reference-frequency branches, `SENSITIVITY_ONLY`.
4. C-03 microwave-index branches, using the Fig.2(b) value once item 1 exists, `SENSITIVITY_ONLY`.
5. Track B single-variable scans with a frozen Track A parent.

Waiting for user authorization: SC-01 cross-section FEM (L2), SC-02 feedline/pad network (L3), SC-03 measured Figure 6/7 comparison.

Auxiliary static learning material is in `docs/learning/`, `notebooks/01_static_mzm.ipynb` and `data/simulated/day01_*`. It does not earn a Track A result grade.
