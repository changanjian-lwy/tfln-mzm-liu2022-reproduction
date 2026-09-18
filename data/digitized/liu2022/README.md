# Liu 2022 raster readout

Source: DOI 10.1109/LPT.2022.3178214, Figure 2(a,c), Figure 3(a inset,b,c). All rows are `digitized_paper_calculation`, never measurements or original FEM arrays.

`curves.csv`: native pixel coordinate, physical frequency/value, extraction uncertainty, panel, series and units. `metadata.json`: native image fingerprints, PDF fingerprint, axis calibration and missing/ambiguous-column counts. Sources stay local. Extracted data carry finite raster precision; many adjacent samples are correlated.

No interpolation is performed in this dataset. Figure 2(c) has missing segments where its 50-ohm dashed reference overlaps the signal. Figure 3(c) omits legend-overlap areas and clipped nulls. Only the 0–50 GHz inset of Figure 3(a) is included; the main panel is not yet extracted. Samples begin above DC. No full nm(f) or complex-Z0 curve is supplied.

Reproduction of extraction requires the user's lawfully obtained PDF and native embedded raster images. On the supplied PDF, PyMuPDF `page.get_images(full=True)` and `doc.extract_image(xref)` yield a 1505x1250 JPEG for Fig.2 (PDF p.2) and a 1250x1404 PNG for Fig.3 (PDF p.3). Save them locally as p2_2.jpeg and p3_0.png, then run `python scripts/digitize_a03.py --images /path/to/local/images`. The script needs NumPy, pandas, Pillow and Matplotlib; PyMuPDF is only needed for the optional PDF extraction step. Different PDF versions must be recalibrated rather than blindly reusing pixels.
