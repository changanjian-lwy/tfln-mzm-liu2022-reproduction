# A01 实际输入与复核记录

**2026-10-04 补写（P1 既有证据审查），不是当时的运行记录。** BOUNDARY 与 RESULTS 在提交 `ecf20af`（2026-09-18）同时进入当前目录；BOUNDARY 写明"declared before execution"，文件修改时间 14:39 早于 RESULTS 的 14:41，但 Git 历史不能独立证明先后。

- 命令：`python scripts/run_a01.py`；`python -m unittest discover -s tests -v`。
- 输入：无外部数据文件；参数写在脚本和 `tests/test_traveling_wave.py` 中（L=6 mm、Vg=1 V、Zg=50 Ω、Z0=50 Ω、nm=ng=2.25、α=0）。
- 代码：`scripts/run_a01.py` 最后修改于 `ecf20af`；`src/tfln_mzm/interaction.py`、`termination.py`、`microwave.py` 最后修改于 `e9a796f`。
- 输出：`dc_results.json`（SHA-256 `ff6143eb1aeac76d`）。
- 环境：见 `docs/environment.md`。
- 复核：2026-10-04 在临时副本中重跑，输出逐字节一致；全部 29 个测试通过，其中 B01–B09 仍通过。

P1 分类：**可引用**，等级 `ANALYTIC_LIMIT_PASS`，范围仅限 L0 解析层。
