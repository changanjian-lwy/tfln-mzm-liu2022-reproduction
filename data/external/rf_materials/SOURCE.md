# RF material constants (Track Q, SC-06)

Evidence class: `external_reference` unless marked `assumption`. Liu 2022 gives no material constants. Track Q treats conductors as perfect (quasi-static, high-frequency limit), so no gold conductivity is used; the internal inductance is only bounded from the paper's own Fig.2(a) loss (Q01 BOUNDARY).

## BCB (benzocyclobutene)

The Dow Chemical Company, "Processing Procedures for CYCLOTENE 3000 Series Dry Etch Resins" (PDF titled "Processing Guidelines for 3000 series of Dry-Etch CYCLOTENE Resins", revised February 2005, 8 pages), Table 1. Retrieved 2026-10-05 from https://wiki.nanofab.ucsb.edu/w/images/7/72/BCB-cyclotene-3000-revA.pdf; kept locally (not in Git; vendor document), SHA-256 `fd75274f5d710e25ed2dbcc3a9d909c1ef31c3dbc3e1ce4ee60104ab83d9d01f`. Transcribed:

| Quantity | Value |
|---|---|
| Dielectric constant | 2.65 – 2.50 at 1–20 GHz |
| Dissipation factor | 0.0008 – 0.002 at 1–20 GHz |

- Liu 2022 (p.856) uses a photosensitive BCB, i.e. the CYCLOTENE 4000 series. Its data sheet (DuPont DS2100, kayakuam.com) returned HTTP 403 and was not read; a search-result summary quotes 2.65 (1 kHz–20 GHz) for it. The same BCB polymer is assumed, so the 3000-series range is used: nominal 2.65, endpoint 2.50.
- Fig.2(b) of Liu 2022 is at 100 GHz, beyond the 20 GHz of the data sheet; using the 1–20 GHz range there is an `assumption`.

## LiNbO3, silica, substrate

- LN: the clamped values of `data/external/dc_materials/SOURCE.md` (ε11 44, ε33 27.9). At 100 GHz the crystal is far above its piezoelectric resonances, so the clamped (constant-strain) values apply (`assumption`). Millimetre-wave measurements of thin-film LN were not accessed; Q02 carries a ×0.95 branch on both components (`assumption`).
- Bonding and PECVD SiO2 (3.8–4.5, nominal 3.9) and substrate (crystalline quartz 4.5, fused silica 3.8): same ranges as `dc_materials/SOURCE.md` (`assumption`).
- Air: 1.
