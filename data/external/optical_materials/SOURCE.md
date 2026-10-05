# Optical material data (refractiveindex.info database)

Evidence class: `external_reference`. Used by Track O (SC-05). None of these numbers come from Liu 2022; they fill gaps the paper leaves open and never override a paper value.

- **Retrieved**: 2026-10-05 (Asia/Singapore) with the GitHub API from https://github.com/polyanskiy/refractiveindex.info-database, branch head `c5c2f188e848453def5970e347399d653df2ffc2` (2026-09-04). The database is public domain (CC0 1.0, stated in each file header). Files are copied byte for byte; only the file names were changed.
- **Formula conventions** (refractiveindex.info): `formula 1`: n² − 1 = C1 + Σ C(2i) λ² / (λ² − C(2i+1)²); `formula 2`: n² − 1 = C1 + Σ C(2i) λ² / (λ² − C(2i+1)). λ in µm. `tabulated nk`: λ (µm), n, k.

| File (database path) | Material | Original source | Used as |
|---|---|---|---|
| `LiNbO3_Zelmon-e.yml` (`main/LiNbO3/nk/Zelmon-e.yml`) | congruent LiNbO3, extraordinary, 21 °C | Zelmon, Small, Jundt, JOSA B 14, 3319 (1997) | main-line LN (n_e) |
| `LiNbO3_Zelmon-o.yml` (`main/LiNbO3/nk/Zelmon-o.yml`) | congruent LiNbO3, ordinary, 21 °C | same | main-line LN (n_o) |
| `MgO-LiNbO3_Zelmon-e.yml` (`other/doped crystals/MgO-LiNbO3/nk/Zelmon-e.yml`) | 5 mol% MgO:LiNbO3, extraordinary | same | sensitivity branch only |
| `MgO-LiNbO3_Zelmon-o.yml` (`other/doped crystals/MgO-LiNbO3/nk/Zelmon-o.yml`) | 5 mol% MgO:LiNbO3, ordinary | same | sensitivity branch only |
| `SiO2_Malitson.yml` (`main/SiO2/nk/Malitson.yml`) | fused silica, 20 °C | Malitson, JOSA 55, 1205 (1965) | bonding and PECVD silica; one substrate branch |
| `SiO2-quartz_Ghosh-o.yml` (`main/SiO2/nk/Ghosh-o.yml`) | crystalline quartz, ordinary | Ghosh, Opt. Commun. 163, 95 (1999) | the other substrate branch |
| `Au_Olmon-ev.yml` (`main/Au/nk/Olmon-ev.yml`) | evaporated gold | Olmon et al., PRB 86, 235147 (2012) | T-rail metal at optical wavelengths |

Caveat recorded by the database itself for the MgO:LiNbO3 files: the original publication appears to have interchanged the ordinary and extraordinary coefficients, and the database swaps them back. Because of that source ambiguity the MgO branch is reported separately and never enters a Track O verdict.

Liu 2022 states only "600-nm-thick X-cut TFLN provided by NanoLN" (p.856); it does not say whether the film is congruent or MgO-doped, nor whether the "quartz" substrate is crystalline or fused. Those choices are declared branches in each Track O boundary.
