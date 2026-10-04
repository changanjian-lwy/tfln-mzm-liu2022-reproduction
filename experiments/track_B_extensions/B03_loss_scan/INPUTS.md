# B03 实际输入与运行记录

边界提交 `30cdcfd`（15:08:21，与 B02–B04 同一次提交）；运行于约 15:08:58（`date` 输出），记录写于 15:10（Asia/Singapore）。

- 命令：`python scripts/run_b02_b04.py --study B03`（复用 `run_b01.inputs` 与 `run_a07.smoothed`）。
- 冻结父实验：A04（`8741206`）；输入 `data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）；G1 对照 B01 `sweep.csv`。
- 输出：本目录 `results.json`、`B03_loss_scan.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：无代码改动。写 RESULTS 后对引用数字逐一核对。
