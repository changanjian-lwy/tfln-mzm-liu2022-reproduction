# A05 convention audit

Parent: A04. Diagnostic only, no fitted parameter promotion and no changes to A03/A04.

Questions: (1) Does the absolute-voltage inset and normalized EO panel imply a consistent low-frequency reference? (2) Can the paper S11 samples satisfy the phase-independent reflection bounds of the current real-Z0 uniform line? (3) Could a mere lossless reference-plane shift explain magnitude residuals?

Fixed inputs: digitized curves, paper L=6mm, Zg=50 ohm, nominal loads. Bracketed input interpolation only, <=1 GHz gaps. No nm/length/loss fit. Calculate implied Vref(f)=|Vavg_fig3a(f)|/10^(M_fig3b(f)/20), on 2–49 GHz supported intersections, and compare its median/spread with ideal DC divider. Propagate source readout uncertainty approximately. This is a cross-panel consistency check, not a measured DC reference.

For real positive Z0 and resistive load, let r=(Z0-Zg)/(Z0+Zg), a=|rhoL|exp(-2alpha L). The full possible |S11| over ALL round-trip phases is bounded by min/max of |(r+a)/(1+r*a)| and |(r-a)/(1-r*a)|. Test the observed magnitude against these bounds. Outside bounds indicates that changing only microwave phase cannot fix the current input model; it does not prove the paper is wrong. Use target ordinate readout uncertainty first, and separately conservative local source/target frequency windows and source ordinate uncertainty. Publish counts, not an adjusted model.

A matched, lossless port extension multiplies S11 by a unit phasor: magnitude cannot change. A mismatched/lossy feedline can change magnitude and requires its own model/data.

No acceptance claim for complete reproduction. Results identify ruled-out explanations and unresolved hypotheses. Do not relabel plot legends, calibrate reference voltages into the baseline, or infer author errors from an incomplete provider.

Secondary diagnostic after reference inference: subtract the inferred reference offset from the A04 EO prediction and report resulting MAE, without evaluating pass/fail. This is explicitly output-calibrated; it is not an independent validation or an adopted correction.
