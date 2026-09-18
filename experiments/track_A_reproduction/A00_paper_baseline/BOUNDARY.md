# A00: paper-baseline preflight

- Parent: initial Git revision 4b08a9d; its static model is auxiliary, not a paper-reproduction result.
- Question: are the declared inputs sufficient to run a numerical Figure 3 reproduction without silently inventing microwave arrays?
- Change: establish a paper-specific input contract and missing-data report. No RF physics run.
- Fixed: Eq.(1)–(4), Figure 3 targets, paper L/Vg/Zg/load values; input metadata in `configs/paper_baseline.json`.
- Success: explicit missing list; no substituted microwave curves; no paper-reproduction pass claim.
- Failure: unknown values silently become constants, units are mixed, sources are absent, or configuration validity is called physics validity.
- Cannot prove: average-voltage correctness, S11 correctness, EO bandwidth, Figure 3 agreement, or any experimental behavior.
- Next experiment: A01 analytical-limit implementation, with explicitly idealized inputs and B01–B09 checks from `docs/BOUNDARIES.md`. A01 is allowed despite missing FEM arrays because its claim is analytical-limit verification only.
