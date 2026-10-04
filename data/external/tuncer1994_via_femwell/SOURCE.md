# Tuncer 1994 CPW reference points (via femwell)

Evidence class: `external_reference`. Used only by A17 (SC-01 stage 1 method gate). These are not Liu 2022 data and never enter a Figure 3 comparison.

- **What**: 14 microwave-index points and 14 attenuation points (dB/cm) for a coplanar waveguide on GaAs, against log10(frequency/GHz).
- **Transcribed from**: femwell, `docs/electronics/examples/RF_CPW_transmission_line_tutorial/RF_CPW_transmission_line_tutorial.py`, arrays `data_nu` (line 1795) and `data_alpha` (line 1815), commit `be2c54763f6215757fd2dac543c9385a82ca5b26` (2025-10-08), https://github.com/HelgeGehring/femwell (GPL-3.0). Only the 28 numbers were copied, not code. A script diff against the tutorial confirmed the transcription is exact.
- **Original source**: E. Tuncer, B.-T. Lee, M. S. Islam, D. P. Neikirk, "Quasi-static conductor loss calculations in transmission lines using a new conformal mapping technique," IEEE Trans. Microw. Theory Tech. 42(9), 1807–1815 (1994), DOI 10.1109/22.310592.
- **Digitized by**: the femwell tutorial authors, not this project. Their digitization error is not stated.
- **Measured or calculated**: the tutorial calls the curve "the measurements from Tuncer1994" and labels it "Tuncer et al. 1992" in its plot. This project has not seen the original figure (paywalled), so whether the points are measurements or the paper's conformal-mapping calculation is **unverified**.
- **Structure parameters** (from the tutorial, also unverified against the paper): signal 7 µm, gap 10 µm, ground 100 µm, metal thickness 0.8 µm, metal conductivity 6×10⁵ S/cm, substrate relative permittivity 13, no dielectric loss.
- **Not copied**: the tutorial's TU Delft simulator curves (`support/*_delft.txt`); the tutorial itself says their Z0 is not expected to match.
