# Source and parameter ledger

Never label generated values as experimental measurements. PDF page numbers below are one-based; journal page numbers are provided separately.

| Category | Label used in project | Interpretation |
|---|---|---|
| Course equation/parameter | `lecture` | Stated in EEK5103 slides; simplified educational model |
| Paper calculation/input | `paper_calculation` / `paper_input` | A calculation or model parameter reported by the authors |
| Paper measurement | `paper_measurement` | A measurement reported by the authors, not this project |
| Digitized curve | `digitized_paper_calculation` or `digitized_paper_measurement` | Samples read from a figure, with extraction uncertainty |
| Assumed parameter | `assumption` | A transparent modeling choice, not recovered author data |
| Our output | `simulation` / `derived_dc_limit` | A result of this project's code or derivation |

## Current ledger

| Quantity | Value | Unit | Category and source | Use |
|---|---:|---|---|---|
| λ0 | 1550 | nm | lecture, Lecture 4 p.22 | Day 1 |
| n_e | 2.2 | 1 | lecture, Lecture 4 p.22 | Day 1 phase index |
| r33 | 30 | pm/V | lecture, Lecture 4 p.22 | Day 1 |
| gap | 10 | µm | lecture, Lecture 4 p.22 | Day 1 |
| length | 20 | mm | lecture, Lecture 4 p.22 | Day 1 |
| Γ | 1 | 1 | assumption consistent with E=V/d lecture example | Day 1 uniform field |
| ideal split and equal loss | 50:50 | — | assumption | Day 1 |
| paper electrode length | 6 | mm | paper_input, Fig.1(a), PDF p.2 / journal p.855 | later RF model |
| paper electrode gap | 3 | µm | paper_input, PDF p.2 / journal p.855 | context only |
| Zg | 50 | Ω | paper_input, PDF p.2 / journal p.855 | source impedance |
| Vg | 1 | V | paper_input, paragraph before Fig.3, p.855 | Fig.3 comparison |
| ZL series | 20, 40, 50, 80, open | Ω | paper_calculation, Fig.3(a,b), p.856 | validation targets |
| n_g | approximately 2.25 | 1 | paper reported optical group index, p.857 | later RF baseline; document any constant-index assumption |
| Vpi | 2.7 | V | paper_measurement, Fig.5, p.856; 1 MHz drive, open terminal | contextual static target only |
| Vpi L | 1.62 | V·cm | authors' reported product, p.856 | contextual comparison |
| measured terminal | 38 | Ω | paper_measurement, Fig.7 discussion, p.857 | separate from Fig.3 40 Ω |
| flatness | <1 over frequencies up to 50 GHz | dB | paper_measurement, Fig.7 discussion, p.857 | not a completed project validation |
| electrical bandwidth | 160 at 6.4 dB | GHz | paper_measurement, Fig.6 discussion, p.856 | not a 3 dB EO bandwidth |

## Missing inputs and extraction policy

The inspected paper includes Eq.(1)–(4), Figure 2 FEM curves and Figure 3 calculated curves. No numerical FEM arrays or field-overlap dataset were found in the downloaded PDF. Z0(f) is not exactly 50 Ω: Figure 2(c) includes oscillations. Do not invent its samples or replace it by 50 Ω without the `assumption` label. Figure 2(a) uses dB/mm; conversion to amplitude attenuation requires α[Np/m] = loss[dB/mm] × 1000 × ln(10)/20.

Figure 3 spans 0–200 GHz. It is calculated; Figures 5–7 contain experimental results and, in Figure 7, modeled extrapolation. No curve samples have been digitized yet. Visual inspection of legends is not a digitized dataset. Later samples must record source panel, load, frequency, ordinate, axis calibration, extraction tool/method, uncertainty and a source-image fingerprint. Avoid presenting sub-pixel precision as accurate data.

## Source inventory

`references/source_inventory.json` records filenames, page counts and SHA-256 fingerprints for all six course lecture PDFs and the downloaded four-page paper. Original PDFs remain outside the repository; neither those PDFs nor course assignments are part of the public project. It must be possible to run the static model without access to the original files.
