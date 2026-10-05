# Track M: measured figures (SC-03)

Authorized 2026-10-05 10:06 (decision D7). Track M reads the paper's **measured** figures (Fig.6, Fig.7) as a separate evidence class, `digitized_paper_measurement`. Nothing here is scored together with the Figure 3 calculation comparison, and nothing here may be used as a Figure 3 input or target ([BOUNDARIES](../../docs/BOUNDARIES.md) §12 item 5). Measured S-parameters are taken at the author's probes and include feedlines and pads; they are not de-embedded line quantities.

Each experiment follows [EXPERIMENT_TEMPLATE.md](../EXPERIMENT_TEMPLATE.md) and commits its BOUNDARY.md before running. Rasters are the native images that `pdfimages -all` extracts from the paper PDF; they stay in the local raster folder (`--images`) and are not in Git.

| Experiment | Grade | Claim boundary |
|---|---|---|
| [M01 Fig.6 digitization](M01_fig6_digitization/RESULTS.md) | PARTIAL_EXTRACTION_COMPLETE | Measured S21 (699 points, 2.4–169 GHz) and S11 (684 points; 12 columns clipped at −36 dB) of the open-terminated device at the probes, and the extracted microwave index (712 points from 0.95 GHz; 35.5% indistinguishable from the red ng = 2.25 line). Calibration matches the author's −6.4 dB line (−6.378 dB) and optical line (2.2502) |
