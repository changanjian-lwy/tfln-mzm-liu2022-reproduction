# M02 实际输入与运行记录

边界提交 `2b80199`（在 M01 运行之前）；脚本提交 `f364ea9`（在 M01 运行之前）；M01 结果与登记 `022a43b`。运行于 2026-10-05 10:23（Asia/Singapore），墙钟 1.0 s。

- 命令：`python scripts/run_m02.py`；事后观察：`python scripts/m02_observation.py --images ../tmp/pdfs/a03`。
- 输入：`data/digitized/liu2022/m01_curves.csv`（SHA-256 `cdbb25e45e5bc5e8`，与 M01 登记值一致）、`data/digitized/liu2022/curves.csv`（`ddf19d3509412656`，A03 Fig.2(a)）；事后观察另读本地栅格 `p4_0.jpeg`（`79198a8026bc874f`）。
- 输出：本目录 `results.json`（`9050d88d47d3a3c2`）、`t3_loss_comparison.csv`（`93ed6ccdedf40a77`）、`m02_checks.png`（`b11f7bd4067a0978`）；事后观察 `observation.json`（`0973c4550e57ee82`，只做静态核对，不重放，与 B05 的事后观察同样处理）。已登记 `experiments/EVIDENCE_MANIFEST.json`。
- 确定性：连续运行两次，三个主输出逐字节相同。
- G1：50 个测试通过；`check_evidence.py` 静态通过（26 个登记实验）；脚本内核对 M01 已登记、哈希一致、M01 判据 2/3/4 对 S21、S11、nm 都通过。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2。
- 实现细节（运行前已写在脚本里，边界没有逐字规定）：T3 插值的不确定度取两侧样本中较大的一个。
- 偏离：无。
