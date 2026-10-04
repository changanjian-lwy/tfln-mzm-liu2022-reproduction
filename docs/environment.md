# Validated local environment

Created 2026-09-18. Project-local `.venv`; no system package changes.

- Python: 3.14.6
- Platform: Darwin 25.6.0, arm64

| Package | Installed version |
|---|---|
| numpy | 2.5.3 |
| scipy | 1.18.1 |
| matplotlib | 3.11.2 |
| pandas | 3.0.6 |
| jupyterlab | 4.6.3 |
| nbformat | 5.11.1 |
| nbclient | 0.11.0 |

`requirements-lock.txt` records installed third-party dependencies. Source PDFs are not runtime dependencies. Scientific libraries and JupyterLab are installed; commercial field solvers are not required.

## Optional FEM layer (added 2026-10-04 for SC-01 / A17)

The cross-section solver `src/tfln_mzm/cpw_fem.py` needs femwell (GPL-3.0, used as a dependency, no code copied), pinned to commit `be2c547` (2025-10-08), plus the packages in `requirements-fem-lock.txt`: scikit-fem 12.0.2, gmsh 4.15.2, pygmsh 7.1.17, meshio 5.3.5, shapely 2.1.2. Install commands are in that file. femwell is installed with `--no-deps` because its declared dependency meshwell pulls in cadquery, vtk and trame, which the solver does not use. Adding these packages left every version in `requirements-lock.txt` unchanged. Without them, `tests/test_cpw_fem.py` is skipped, and `check_evidence.py --replay` skips any experiment whose replay entry declares `requires: ["femwell"]`. A17 itself is registered without a replay entry: its eigen-solves differ between runs at the 1e-4 relative level, so its outputs are hash-checked only.
