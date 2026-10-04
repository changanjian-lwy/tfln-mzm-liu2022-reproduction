# A14 实际输入与运行记录

边界提交 `44bbc9b`（15:11:59）；运行于 15:12:32（`date` 输出），记录写于 15:13（Asia/Singapore）。

- 命令：`python scripts/run_a14.py`。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）、`a06_curves.csv`；A04 与 A13 的 `results.json`（存档 MAE 与通过状态）。
- 输出：本目录 `results.json`、`corner_envelopes.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：无。一次运行即得到本结果。
