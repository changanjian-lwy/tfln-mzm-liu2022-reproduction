# A07 实际输入与运行记录

边界提交 `42f09e0`；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 15:40（Asia/Singapore）。

- 命令：`python scripts/run_a07.py`。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）；A04 `results.json`（只读取存档 MAE 作 G1 对照）；Fig.2(b) 星标值 202.5997 GHz 取自 A06 `results.json`，作为常数写在脚本中。
- 输出：本目录 `results.json`、`bandwidth_definitions.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：第一次运行前发现脚本把模型频率以 Hz 传入、而平滑窗口与缺口判据按 GHz 写，改为全程 GHz 后才产生本结果；判据、窗口、阈值都没有改。没有任何以 Hz 错误单位得到的结果被保存或用于判断。
