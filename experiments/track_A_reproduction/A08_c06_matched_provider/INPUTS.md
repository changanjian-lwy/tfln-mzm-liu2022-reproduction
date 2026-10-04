# A08 实际输入与运行记录

边界提交 `f714833`；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 16:05（Asia/Singapore）。

- 命令：`python scripts/run_a08.py`。
- 输入：`data/digitized/liu2022/curves.csv`（Fig.2(a) 损耗，SHA-256 `ddf19d3509412656`）、`a06_curves.csv`（Fig.2(b) 两条曲线）、A06 `results.json`（星标与虚线读数）。
- 假设（分支 C-06b）：Z0=Zg=ZL=50 Ω；损耗不随 BCB 变化；nm 无色散；Fig.2(a) 覆盖区外用加权最小二乘 α=a√f+bf（a=0.03088 dB/mm/√GHz，b=0.001906 dB/mm/GHz，残差 0.012 dB/mm）。
- 输出：本目录 `results.json`、`matched_provider.png`；登记在 `experiments/EVIDENCE_MANIFEST.json`。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：第一次运行在写 JSON 时因 numpy 布尔类型报错退出，没有产生任何结果文件；把该项包进 `bool()` 后重跑。判据、输入、假设都没有改。另在提交前单独核实了 RESULTS 限制 3 的说法：nm 偏 ±0.0005 时设计点响应变化约 1e-4 dB。
