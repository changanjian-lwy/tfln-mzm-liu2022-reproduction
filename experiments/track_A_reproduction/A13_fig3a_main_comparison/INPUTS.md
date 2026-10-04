# A13 实际输入与运行记录

边界提交 `2e94547`（15:05:57）；运行于 15:06:21（`date` 输出），记录写于 15:07（Asia/Singapore）。

- 命令：`python scripts/run_a13.py`。
- 输入：`data/digitized/liu2022/a06_curves.csv`（Fig.3(a) 主图读数，哈希写在 `results.json`）、`curves.csv`（Fig.2 输入与插图，SHA-256 `ddf19d3509412656`）；A04 `results.json`（存档 MAE，G1 对照）。
- 输出：本目录 `results.json`、`fig3a_main_comparison.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：无。一次运行即得到本结果。
