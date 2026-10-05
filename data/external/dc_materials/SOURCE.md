# Low-frequency material constants (Track O, O02)

Evidence class: `external_reference` unless marked `assumption`. None of these come from Liu 2022, which gives no material constants. The O02 boundary uses only the values below.

## LiNbO3 (one provider)

Crystal Technology, Inc., "Lithium Niobate" optical-crystal data sheet (PDF titled "DC12_LN OptCrys_8/99", 2 pages, 1999), formerly at www.crystaltechnology.com. Retrieved 2026-10-05 from a third-party mirror; kept locally (not in Git; vendor document), SHA-256 `a83393d0ef69c9a92133668b5a9b1a47fe6a83f80300156ece57cfe6338e384d`. Transcribed values (25 °C):

| Quantity | Unclamped ("n < 500 kHz") | Clamped ("n > 10 MHz") |
|---|---|---|
| ε11 (crystal X, Y) | 85 | 44 |
| ε33 (crystal Z) | 28.7 | 27.9 |
| r33 at 633 nm (pm/V) | 33 | 31 |

- The sheet gives r33 only at 633 nm. The value at 1550 nm is not in hand; normal dispersion makes it lower, which would raise VπL. O02 does not cover that direction and says so.
- Ref [10] of Liu 2022 (Liu et al., Chin. Opt. Lett. 19, 060016, 2021, p.1) uses r33 = 31 pm/V, the same as the clamped value here.
- Weis & Gaylord, Appl. Phys. A 37, 191 (1985) and Jazbinšek & Zgonik, Appl. Phys. B 74, 407 (2002) are the usual primary compilations; this project has not accessed them.
- The 1 MHz drive used for Vπ in Liu 2022 lies between the sheet's two regimes, so clamped and unclamped are both declared branches. Each branch pairs its own ε and r33.

## Silica and substrate (assumption ranges)

No source in hand pins these for the Liu device; the ranges below span the values met while searching and are labelled `assumption`:

- Bonding and PECVD SiO2: 3.8–4.5, nominal 3.9. Pointers: fused-quartz vendor sheets quote 3.8; MIT 6.777 material-property page lists 5 for gap-filling PECVD SiO2 (Schwartz et al., J. Electrochem. Soc.), outside the range and not adopted; reports of as-deposited PECVD SiO2 give 3.84–4.2.
- Substrate: crystalline quartz 4.5 (reported values 4.3–4.6, anisotropic) or fused silica 3.8, tied to the O01 substrate branch.
- Air: 1.
