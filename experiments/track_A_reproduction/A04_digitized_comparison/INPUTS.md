# A04 实际输入与复核记录

**2026-10-04 补写（P1 既有证据审查），不是当时的运行记录。** BOUNDARY 与 RESULTS 在提交 `8741206`（2026-09-18）同时进入 Git。BOUNDARY 写有事前判据（≥100 个评分点、≥90% 落在读数容差内、无未解析的非有限预测），但 Git 历史不能独立证明它先于结果写定。

- 命令：`python scripts/run_a04.py`。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`，同时写在 `results.json` 的 `source_sha256`）。
- 代码：`scripts/run_a04.py` 最后修改于 `8741206`。
- 输出：`comparison.csv`（`63f501a927d2df03`）、`comparison.png`（`c8090b704a22420f`）、`results.json`（`dce5c94dfcfd3a28`）；`results_before_axis_tick_QA.json`（`4947ba0fb05126b2`）是 A03 轴刻度 QA 之前的历史报告，不能重新生成，只做静态核对。
- 偏离（当时已披露）：初次对照后 A03 修正了轴刻度掩码并重新生成；物理参数与判据未改，评分结果未变。
- 环境：见 `docs/environment.md`。
- 复核：2026-10-04 在临时副本中重跑，三个可再生输出逐字节一致。

P1 分类：**可引用**，等级 `PARTIAL_COMPARISON_FAIL`，层级 L1。不能引用为 Figure 3 复现。
