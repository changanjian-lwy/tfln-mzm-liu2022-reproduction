# A11 实际输入与运行记录

边界提交 `775d1e9`；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 17:45（Asia/Singapore）。

- 命令：`python scripts/run_a11.py`（从 `run_a08.py` 导入损耗输入侧拟合，从 `run_a09.py` 导入读数与负载辅助函数）。
- 输入：`data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`）；A05 `results.json`（隐含偏移，只作比较）；A04 `results.json`（存档 MAE，G1 对照）。
- 假设：2.22 GHz 以下 Z0 保持 Fig.2(c) 第一个读数（46.29 Ω），损耗用 α=a√f+bf 拟合（A08 同一拟合）。
- 输出：本目录 `results.json`、`reference_frequency.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：第一次运行后加入一项只报告的保持值敏感性（±2 倍中位 Z0 读数误差），重跑。候选、判据、容差都没有改，主结果两次相同。
