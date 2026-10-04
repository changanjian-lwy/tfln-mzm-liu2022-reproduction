# A00 实际输入与复核记录

**2026-10-04 补写（P1 既有证据审查），不是当时的运行记录。** BOUNDARY 与 RESULTS 都在提交 `ecf20af`（2026-09-18）首次进入当前目录（更早的配置在 `884a480`）；两者同一提交，事前声明无法由 Git 历史独立证明。

- 命令：没有独立脚本；`results.json` 等于 `tfln_mzm.contracts.preflight(configs/paper_baseline.json)` 的输出，由 `tests/test_contracts.py` 覆盖。
- 输入：`configs/paper_baseline.json`（SHA-256 `b2553ee542fc8b1c`）。
- 输出：`results.json`（`f80b2da2d37b0e61`）。
- 环境：见 `docs/environment.md`（Python 3.14.6、numpy 2.5.3）。
- 复核：2026-10-04 在临时副本中重算 preflight，结果与存档一致；全部 29 个测试通过。

P1 分类：**可引用**，作为输入缺失清单的历史快照。注意配置中 Z0、α 仍标为 `unknown`，因为 A00 早于 A03 的数字化；这是当时的正确状态，不需要改写。
