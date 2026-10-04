# B06 实际输入与运行记录

边界提交 `7e6e136`（00:37:05），脚本提交 `f7f8909`（00:42:16）；运行于 2026-10-05 00:42:16–00:42:27（Asia/Singapore，取自 `date` 与 git）。

- 运行前 G1：`python -m unittest discover -s tests` 50 项通过；`python scripts/check_evidence.py --replay` 通过（含刚登记的 B05，五个输出逐字节相同；A17、A17b 照旧只做静态核对）。
- 命令：`python scripts/run_b06.py`（墙钟 10.7 s，峰值内存约 0.39 GB；预算 600 s）。G2 不过，脚本按设计只写 `results.json`（等级 `NUMERICAL_FAIL` 与判据），不写前沿、BO 和图。
- 失败诊断（看到失败以后才写，只打印、不登记）：`python scripts/b06_g2_diagnostics.py`。
- 输入：`data/digitized/liu2022/curves.csv` `ddf19d3509412656`；`experiments/track_B_extensions/B01_termination_sweep/sweep.csv` `90368de2db7d1a35`；代码复用 `scripts/run_b05.py`（`35c63b2`）。
- 输出：`results.json`，见 EVIDENCE_MANIFEST.json 的 `B06_flat_covered_design_bo` 条目。
- 环境：同 B05；没有新增依赖。
- 可重放性：在 scratch 副本中重跑，`results.json` 逐字节相同。
- 偏离：无。
