# Q01 实际输入与运行记录

边界提交 `3d39437`（16:42:19）；运行前修订、模块、测试与脚本 `94ed1a7`（16:57:54，提交时没有用有限元解过 Liu 的任何截面）。登记的运行：2026-10-05 16:57:59–16:58:48（Asia/Singapore），20 次求解，墙钟约 49 s，单次最长 2.4 s，峰值内存约 1.0 GB，最大网格 99 704 单元；没有触发 `BUDGET.json` 的上限。16:59–17:00 第二次运行：`results.json` 与 `q01_inductance_z0.png` 逐字节相同，只有 `run_log.json`（计时）不同；仓库里的是第二次的文件。

- G1 的非 FEM 部分（脚本提交前）：74 个测试通过（67 个旧测试 + 7 个新测试）；`check_evidence.py` 静态通过（31 个登记实验）。
- 命令：`PYTHONPATH=src .venv/bin/python scripts/run_q01.py`（退出码 0）。
- 输入：`data/digitized/liu2022/curves.csv`（A03，`ddf19d3509412656`：Fig.2(a) 损耗与 Fig.2(c) Z0）；nm = 2.2500 ± 0.0005 写在脚本里，取自 A06 星标（`data/digitized/liu2022/a06_curves.csv` `27c3c45e90d63dcc` 与 A06 RESULTS）；`BUDGET.json`（`7fb2c68a7cb8885b`）。脚本不读任何 PDF。
- 输出：`results.json`（`ea878e15ceee8717`）、`q01_inductance_z0.png`（`0ef7cd9e53b0db3d`）、`run_log.json`（只登记静态哈希）；已登记 `experiments/EVIDENCE_MANIFEST.json`，可重放（需要可选 FEM 依赖）。
- 网格：接受 m = 2。基准 7 728 单元；地宽 150 µm 为 14 890（m = 2）、52 907（m = 1）、25 218（d = 2）；地宽 1000 µm 为 62 574。
- 环境：Python 3.14.6、numpy 2.5.3、scipy 1.18.1、pandas 3.0.6、scikit-fem 12.0.2、gmsh 4.15.2、pygmsh 7.1.17、meshio 5.3.5、shapely 2.1.2、matplotlib 3.11.2、femwell `be2c547`（只用其网格函数）；没有新装依赖。

## 偏离
- 无。运行前的方法修改（外边界改为自然边界、网格尺寸写明）和非 Liu 截面的干跑都在 BOUNDARY 第 5 节披露，提交在运行之前。
