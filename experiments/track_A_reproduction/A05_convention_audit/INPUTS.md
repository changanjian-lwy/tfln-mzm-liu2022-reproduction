# A05 实际输入与复核记录

**2026-10-04 补写（P1 既有证据审查），不是当时的运行记录。** BOUNDARY 与 RESULTS 在提交 `be41c24`（2026-09-18）同时进入 Git，事前声明无法由 Git 历史独立证明。

- 命令：`python scripts/run_a05.py`。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）；A04 的 `comparison.csv`（`63f501a927d2df03`），用于第二部分输出校准诊断。A05 因此依赖 A04 的存档输出。
- 代码：`scripts/run_a05.py` 与 `tests/test_audit.py` 最后修改于 `be41c24`。
- 输出：`phase_bounds.csv`（`c3756c35b61068b9`）、`results.json`（`8816ec296e9adb19`）、`normalization_diagnostic.json`（`097d8f60b4d38730`）。
- 环境：见 `docs/environment.md`。
- 复核：2026-10-04 在临时副本中重跑，三个输出逐字节一致。

P1 分类：**可引用**，等级 `DIAGNOSTIC_ONLY`。它指认的候选原因已登记为冲突分支 C-01、C-02、C-05；输出校准后的 MAE 只是诊断，不是被采用的修正。
