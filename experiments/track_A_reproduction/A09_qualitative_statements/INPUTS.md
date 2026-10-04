# A09 实际输入与运行记录

边界提交 `c4445e9`；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 16:30（Asia/Singapore）。

- 命令：`python scripts/run_a09.py`。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）；A04 `results.json`（只读取存档 MAE 作 G1 对照）。
- 输出：本目录 `results.json`、`qualitative_checks.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：第一次运行时，代码把 S2 的"无法判定"（模型在 10 GHz 不受支持）当成了"不通过"，与 BOUNDARY 第 1 节"少于规定点数记为无法判定"不符。改为返回"无法判定"，并把判定字段改为逐条文字结论后重跑。阈值、窗口、写法都没有改；读数一侧的结果两次完全相同。
