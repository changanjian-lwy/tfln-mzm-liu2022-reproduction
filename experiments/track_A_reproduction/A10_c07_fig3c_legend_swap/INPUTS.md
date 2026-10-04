# A10 实际输入与运行记录

边界提交 `fc54cde`；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 16:55（Asia/Singapore）。

- 命令：`python scripts/run_a10.py`（从 `scripts/run_a05.py` 导入 `bounds`，不重新推导反射界）。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）；A05 `results.json` 与 A04 `results.json`（只读取存档数字作 G1 对照）。
- 输出：本目录 `results.json`、`legend_swap.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：无。一次运行即得到本结果。
