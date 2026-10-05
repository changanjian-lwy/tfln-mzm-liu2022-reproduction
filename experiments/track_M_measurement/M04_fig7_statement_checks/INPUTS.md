# M04 实际输入与运行记录

边界提交 `d432c27`（在 M03 运行前）；脚本提交 `c8f1b39`（在 M03 运行前）；运行前修改与 M03 结果一起在 `cfa76cd` 提交。运行于 2026-10-05 11:32（Asia/Singapore），墙钟 0.5 s。

- 命令：`python scripts/run_m04.py`；事后观察：`python scripts/m04_observation.py --images ../tmp/pdfs/a03`。
- 输入：`data/digitized/liu2022/m03_curves.csv`（SHA-256 `41d8577b9fd0fa05`，与 M03 登记值一致）、`m03_metadata.json`（读截断列表）、M03 `results.json`（读判据）；事后观察另读本地栅格 `p4_1.png`（`6efa75979dca4e92`）。
- 输出：本目录 `results.json`（`906f6bff4a7b0f01`）、`m04_checks.png`（`4220c0573139709a`）；事后观察 `observation.json`（`a52999729349c29e`）与 `OBSERVATION_posthoc.md` 只做静态核对，不重放。已登记 `experiments/EVIDENCE_MANIFEST.json`。
- 确定性：连续运行两次，两个主输出逐字节相同。
- G1：50 个测试通过；`check_evidence.py` 静态通过（28 个登记实验）；脚本内核对 M03 已登记、哈希一致、M03 判据 2/3 通过、38 Ω 的 EO 与 S11 判据 4/5 通过。
- 运行前修改：脚本提交 `c8f1b39` 之后、首次运行之前（已看过 M03 结果）只改了一处——"并报不判"的其他负载只用 M03 判为可用的系列，因为 M03 判据 4 把 EO 83 Ω 标为不可用，并规定下游不得使用。U1a、U1b、U2 的规则与对象（只用 38 Ω）没有改。截断列只计下框（低于 −35 dB）这一条，在 `c8f1b39` 提交之前已写进脚本，与边界"截断点（低于 −35 dB）算作低于"一致。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 偏离：无。
