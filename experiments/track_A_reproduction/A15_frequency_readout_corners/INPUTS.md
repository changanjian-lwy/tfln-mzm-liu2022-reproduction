# A15 实际输入与运行记录

边界提交 `48c6377`（15:16:29）；运行于 15:17:01（`date` 输出），记录写于 15:17（Asia/Singapore）。

- 命令：`python scripts/run_a15.py`。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）、`a06_curves.csv`；A14 `results.json`（G1 对照与原判定）。
- 输出：本目录 `results.json`、`frequency_corners.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：第一次运行输出了 numpy 的 "All-NaN slice" 警告（来自不参与评分、没有任何受支持角点的频点）；加了只针对该警告的屏蔽后重跑，`results.json` 逐字节相同（SHA-256 前 16 位 `16b3b336d4a0d8f3`）。
