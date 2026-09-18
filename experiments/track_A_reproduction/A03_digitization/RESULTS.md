# A03 extraction results

Status: PARTIAL_EXTRACTION_COMPLETE; claim scope is raster digitization, not model validation.

Extracted input curves: Figure 2(a) this-work loss, Figure 2(c) Z0. Extracted independent target curves: Figure 3(a) five-load inset, Figure 3(b) five-load EO response, Figure 3(c) four-load S11. Native raster image and PDF fingerprints, axis calibrations, per-column source pixels and extraction estimates are saved with the dataset.

Native rasters and local extraction overlays were inspected. The first preview exposed a black legend-text descender contaminating Figure 3(a) near 29 GHz; the black-only legend mask was expanded to y>184 and data regenerated. No physical parameter or target ordinate was adjusted to match the model. All masked/multiple-cluster columns are omitted, not interpolated. Fig.2(c)'s longest gap is about 6.67 GHz; later models must declare how such gaps are treated. A stricter red/blue/black intensity threshold produces median common-point differences of zero; this tests color localization only and does not validate axis calibration or recover hidden samples.

Known limitations: no complete 0–200 GHz Fig.3(a), no nm(f), no complex Z0, no exact DC/1 MHz data, clipped deep S11 nulls and discarded legend regions. Extraction errors are estimated from 2 pixels plus stroke thickness, not statistical confidence intervals. Local overlays are excluded from Git because they contain the original paper figure.

Next: freeze the dataset and compare a declared model only within supported input intervals, excluding long Z0 gaps. Low-frequency normalization and constant nm require explicit assumptions. Quantitative acceptance criteria must be declared before evaluating candidate outputs.

Final overlay QA also identified bottom-axis tick fragments in the black inset series at about 30 and 40 GHz. Black pixels with y>=220 are now excluded (the visible black trace is well above that row). A04 pre-QA metrics were retained separately; the final comparison was regenerated with identical physical parameters and unchanged acceptance. Scored metrics were unaffected because those isolated tick samples had already failed target-gap eligibility.
