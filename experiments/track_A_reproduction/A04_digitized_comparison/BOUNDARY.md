# A04: first independent numerical comparison

Parent: A01 equations and frozen A03 readout. A02 illustrative loss is replaced by Figure 2(a) readout; constant Z0 is replaced by Figure 2(c) readout. This is an explicitly declared two-input provider replacement, not a one-variable sensitivity study. No parameter fit.

Fixed: L=6 mm, Vg=1 V, Zg=50 ohm; targets 20/40/50/80 ohm/open. Assumption: nm=ng=2.25 independent of frequency, following nominal matching but not recovering unavailable dispersion. Use the previously chosen f0=1 MHz: alpha(f0)=0 and Z0(f0)=50 ohm as a separate low-frequency approximation, so Vref=Vg*ZL/(Zg+ZL). This approximation is not digitized and is not propagated into the higher-frequency provider.

Evaluate only at observed target frequencies for which both source curves bracket the frequency with gaps <=1 GHz. Linear interpolation only; no extrapolation; long gaps rejected. This excludes some frequency intervals and cannot produce a complete bandwidth result. No smoothing or phase adjustment against targets.

Predeclared diagnostic acceptance per target series: >=100 points and >=90% of retained samples within the target's extraction allowance, defined as uncertainty_value + abs(local target slope)*uncertainty_frequency_GHz. Slope uses neighboring retained target samples and is an approximate horizontal-error propagation. Exclude target samples adjacent to target gaps >1 GHz from this scoring. This is a strict diagnostic of agreement within estimated target readout precision; it excludes unknown input/interpolation and model errors, so failure does not prove the paper is wrong. Report MAE, RMSE, maximum error and fraction within allowance for every series, even if failed. Do not modify acceptance or parameters after inspection.

S11 matched-load null singularities: nonfinite predicted dB values must be reported as unresolved, not clipped or silently dropped. Missing input coverage must be counted separately. Figure 3(a) validation is inset-only. No measured data are used. A result cannot claim full Figure 3 replication unless all panels and missing physical inputs are independently resolved.
