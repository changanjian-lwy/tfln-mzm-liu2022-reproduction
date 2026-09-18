# Experiment index

Read `docs/BOUNDARIES.md` before running a new experiment. Every experiment owns a BOUNDARY.md and RESULTS.md; preserve both successful and failed conclusions. A new physical assumption or provider gets a new experiment rather than overwriting a completed baseline.

| Track | Experiment | Status | Claim boundary |
|---|---|---|---|
| A: paper reproduction | [A00 baseline](track_A_reproduction/A00_paper_baseline/RESULTS.md) | Missing input arrays explicitly recorded | Preflight only |
| A: paper reproduction | [A01 analytical limits](track_A_reproduction/A01_analytic_limits/RESULTS.md) | ANALYTIC_LIMIT_PASS | Eq.(1)–(4) / S11 analytical checks, not Figure 3 agreement |
| A: paper reproduction | A02 frequency-response and bandwidth extraction | Planned | First downward crossing, censored bandwidth, convergence |
| A: paper reproduction | Figure 2 inputs / Figure 3 digitized comparison | Planned | Uncertainty-labelled reproduction |
| B: extensions | Length, RF loss and mismatch scans | Planned | Separate from paper baseline |

Auxiliary static learning material is in `docs/learning/`, `notebooks/01_static_mzm.ipynb` and `data/simulated/day01_*`. It does not earn a Track A result grade.
