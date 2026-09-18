# Day 1 validation record

Date: 2026-09-18. Status: passed.

- Seven physical/unit checks passed with `python -m unittest discover -s tests -v`.
- The lecture benchmark agrees with manual substitution: Vpi = 2.426120711 V, |delta n| = 3.875e-5.
- Independent complex-field interference agrees with the intensity expression and conserves total output power.
- Half-wave switching, periodicity, quadrature slope, length/gap/overlap scaling, push-pull convention and DC divider limits pass.
- The notebook executed all code cells successfully using the project environment.
- The deterministic script generated the figure, 1001-point static table, DC-limit table and JSON summary.
- Figure checked visually for labels, units, layout and simulation identification.
- Installed package dependency check reports no broken requirements.

Limits: these checks validate the implementation of the ideal static equations, not the physical accuracy of a fabricated device. No RF/EO bandwidth, Figure 3 error metric, or new experimental result is established by this record. The notebook is a teaching artifact; tests are the independent checks of the reusable functions.
