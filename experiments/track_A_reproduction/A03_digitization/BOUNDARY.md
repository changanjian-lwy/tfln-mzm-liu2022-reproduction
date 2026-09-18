# A03: independent raster curve extraction

Parent: A02. Scope: extract published calculated curves, not validate the model.

Input: locally supplied four-page DOI 10.1109/LPT.2022.3178214 PDF. Embedded Figure 2 and Figure 3 raster images are extracted without rescaling; source SHA-256 fingerprints recorded. No source images are redistributed. No model results enter extraction.

Before extraction, lock linear-axis calibration in native image pixels: Fig2a x=120..687 -> 0..200 GHz, y=526..106 -> 0..1.5 dB/mm; Fig2c x=495..1035 -> 0..200 GHz, y=1150..732 -> 0..100 ohm; Fig3b x=195..1111 -> 0..200 GHz, y=858..522 -> -12..6 dB; Fig3c x=195..1111 -> 0..200 GHz, y=1320..974 -> -50..0 dB. Fig3a inset x=650..1107 -> 0..50 GHz, y=258..24 -> 0.2..1 V. Coordinates selected by visual inspection of native pixels.

Method: color discrimination and per-column connected row clusters, excluding axes and legends. Retain only a single cluster per column; ambiguous/missing columns are omitted and counted, not filled. Fig2c's dashed 50-ohm reference cannot be reliably separated where it overlaps the curve: exclude y=939..947 and do not bridge resulting gaps. Fig3c legend excluded at x>=740,y>=1180. Fig3a inset legend excluded above y=172 for x>=840 (black trace uses y=184 after visual QA found a legend descender at y≈173). These masks limit data coverage and must stay recorded.

Uncertainty: raster localization estimate ±2 pixels on each axis, enlarged vertically by half detected stroke thickness. This is an extraction estimate, not confidence bounds on author simulations, and does not include an unknown interpolation/model error. Compare separately with a stricter color threshold; median ordinate differences and valid sample counts are recorded. Visually inspect the extracted curves and local overlays. Low/high edge pixels excluded by 4 px. No point below available frequency support or exact DC sample is invented.

Success: source fingerprints, calibrations, masks, provenance and finite samples retained; extraction overlays visually acceptable; missing regions explicit. Failure: selecting a legend/reference line as signal, silently filling gaps, or claiming that raster readout is original FEM data. Never fit the model while extracting targets.

Not proven: complete Figure 3(a) 0–200 GHz coverage, full complex Z0, nm(f), low-frequency extrapolation, or quantitative model agreement. Only inset Figure 3(a) is extracted in this experiment; the full panel remains future work.

QA mask addendum: black inset samples y>=220 are excluded as bottom-axis tick artifacts; this is source-image QA, not a fitted signal constraint.
