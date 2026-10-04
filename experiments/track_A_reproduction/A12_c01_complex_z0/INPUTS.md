# A12 实际输入与运行记录

边界提交 `224ef3a`；模块扩展与新测试 `a2890b1`（实验运行前提交）；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 18:20（Asia/Singapore）。

- 命令：`python scripts/run_a12.py`；模块关卡 `python -m unittest discover -s tests` 与 `python scripts/check_evidence.py --replay`。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）；A04 `results.json`（存档 MAE，G1 对照）。
- 代码改动：`src/tfln_mzm/termination.py` 新增复数 Z0 路径；`tests/test_complex_z0.py` 新增 5 个测试。新测试在旧模块上运行时 4 个报错（旧模块会把复数 Z0 截成实数），在新模块上全部通过。
- 输出：本目录 `results.json`、`complex_z0.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、scipy 1.18.1、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：模块扩展时，"所有实验目录都必须登记"的测试会让"只交了边界、尚未运行"的实验失败；改为"有 RESULTS.md 的实验必须登记、每个目录必须有 BOUNDARY.md"。这是证据规则的细化，不涉及本实验判据。
