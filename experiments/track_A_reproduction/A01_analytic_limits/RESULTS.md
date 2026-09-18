# A01 results

Grade: `ANALYTIC_LIMIT_PASS`. Scope: uniform effective-line analytical implementation, Eq.(1)–(4) magnitude interpretation plus derived S11. No Figure 3 numerical agreement is claimed.

Execution completed normally. The full suite passed 20 tests: 10 new traveling-wave/unit-domain tests, three preflight tests, seven auxiliary static tests. B01–B09 all passed their predeclared tolerances. The B09 oracle directly integrates the printed Eq.(1) in a separate coordinate; the S11 oracle uses ABCD matrices rather than the reflection helper.

Parameters and variations are declared in BOUNDARY.md and tests/test_traveling_wave.py. `python scripts/run_a01.py` regenerates the DC table and actual inputs in dc_results.json. `python -m unittest discover -s tests -v` repeats validation.

No failure boundary encountered within the selected test domain. Missing original FEM arrays remain missing in A00; the present idealized tests do not resolve that evidence gap. Real positive Z0 and resistive nonnegative loads (including open) define the current solver domain. Complex characteristic impedance and reactive loads are not implemented. B10–B12 bandwidth logic remains unimplemented.

Next: A02 frequency response and first-crossing/censored bandwidth extraction under a separately frozen, explicitly assumed provider. After that, digitized paper inputs and outputs may be compared. Do not adjust the paper loads, source voltage or source impedance to force a fit.
