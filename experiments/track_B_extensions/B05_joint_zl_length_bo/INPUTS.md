# B05 实际输入与运行记录

边界提交 `46ade0c`（00:19），模块与脚本提交 `35c63b2`（00:24）；运行于 2026-10-05 00:24:29–00:28:18（Asia/Singapore，取自 `date`）。

- 运行前 G1：`python -m unittest discover -s tests` 50 项通过（原 44 项 + `tests/test_surrogate.py` 6 项；交接里写的"46 项"有误，`30fe98b` 时是 44 项）；`python scripts/check_evidence.py --replay` 通过（A17、A17b 照旧只做静态核对）。
- 命令：`python scripts/run_b05.py`（墙钟 229 s，峰值内存约 1.5 GB；预算 900 s）。
- 事后观察（运行后才写，不是判据）：`python scripts/b05_observation_peaking.py`，只读 `front.csv`，用同一提供者重算曲线。
- 输入：`data/digitized/liu2022/curves.csv` `ddf19d3509412656`；`experiments/track_B_extensions/B01_termination_sweep/sweep.csv` `90368de2db7d1a35`；`experiments/track_B_extensions/B02_length_scan/results.json` `c5a45f3d8eb1b1ce`（与 EVIDENCE_MANIFEST 一致）。
- 输出：见 EVIDENCE_MANIFEST.json 的 `B05_joint_zl_length_bo` 条目（`results.json`、`front.csv`、`bo_runs.csv`、`b05_design_space.png`、`b05_front_bo.png`；事后观察 `observation_peaking.json`、`observation_response.png`）。
- 环境：Python 3.14.6、numpy 2.5.3、scipy 1.18.1、pandas 3.0.6、matplotlib 3.11.2；没有新增依赖。
- 可重放性：在 scratch 副本中重跑 `run_b05.py`，五个输出与工作区逐字节相同；事后观察脚本重跑后两个输出也相同。登记时 `run_b05.py` 作为重放命令；两个观察输出登记为 static_only（只核对哈希；单个实验只能登记一条重放命令）。
- 偏离：(1) 事后观察脚本是看到结果后加的，输出单独成文件，不进入判据；(2) `b05_front_bo.png` 没有标出 B=50/75 GHz 的最优落在 L=30 mm 边界（截尾），在 RESULTS 第 0 节第 2 条说明，没有重画；(3) 没有其他偏离。
