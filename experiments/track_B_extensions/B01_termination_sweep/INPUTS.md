# B01 实际输入与运行记录

边界提交 `fff5290`；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 17:15（Asia/Singapore）。

- 命令：`python scripts/run_b01.py`（从 `scripts/run_a09.py` 导入已验证的 A04 复算函数）。
- 冻结父实验：A04（`8741206`）；输入 `data/digitized/liu2022/curves.csv`（SHA-256 `ddf19d3509412656`），A04 `results.json`（只读取存档 MAE）。
- 输出：本目录 `results.json`、`sweep.csv`、`termination_sweep.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：第一次运行后，在脚本中加入一项只报告的观察量（同一起伏定义在作者 Fig.3(b) 读数上的值），重跑。判据、扫描、输入都没有改，模型结果两次相同。写 RESULTS 后对引用的数字逐一核对，改正了两处（S11 区间 25–67 Ω、驱动电压代价 0.5–1.3 dB）。
