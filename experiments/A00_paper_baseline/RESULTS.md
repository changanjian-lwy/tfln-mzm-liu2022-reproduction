# A00 results

Configuration preflight completed normally. Its missing-data guards passed.

The Figure 3 numerical reproduction remains `BLOCKED_BY_MISSING_DATA` for the undeclared Z0(f), attenuation(f), and microwave-index(f) providers. This does not block analytical-limit tests with separately declared idealized inputs. The exact missing list is in `results.json`.

Actual inputs: `configs/paper_baseline.json`. No propagation simulation or fitting was run. First blocking boundary: unavailable microwave input providers. The reference frequency is an explicit project assumption. Ten tests pass: three configuration-guard tests and seven auxiliary static tests; none is a passed RF boundary test.

Next: A01 implements analytical limits before Figure 3 fitting. Do not change the paper L, Vg, Zg or Figure 3 load series to remove the missing-data status. Each data provider must be supplied with its own evidence and applicable frequency range.
