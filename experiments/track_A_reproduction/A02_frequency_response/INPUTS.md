# A02 实际输入与复核记录

**2026-10-04 补写（P1 既有证据审查），不是当时的运行记录。** BOUNDARY 与 RESULTS 在提交 `35d831b`（2026-09-18）同时进入 Git，事前声明无法由 Git 历史独立证明。

- 命令：`python scripts/run_a02.py`。
- 输入：无外部数据文件；假设提供者写在脚本中（Z0=50 Ω、nm=ng=2.25、损耗 0.05·√(f/GHz) dB/mm、f0=1 MHz、网格 2001/4001 点）。
- 代码：`scripts/run_a02.py` 与 `src/tfln_mzm/response.py` 最后修改于 `35d831b`。
- 输出：`curves.csv`（SHA-256 `f4a0f2f2613e712d`）、`response.png`（`90ea442f1cf1ed25`）、`results.json`（`4145ea6076ae800c`）。
- 环境：见 `docs/environment.md`。
- 复核：2026-10-04 在临时副本中重跑，三个输出逐字节一致（含 PNG）。

P1 分类：**可引用**带宽提取逻辑（B10–B12）；物理数字只是假设损耗律下的机理演示，等级 `QUALITATIVE_ONLY`，不能当成论文器件带宽。
