# A17b 实际输入与运行记录

边界提交 `1eead28`（22:36:08）；运行前修订 `f962b12`（22:46:57，G1 未通过）；模块、测试与（未运行的）脚本 `5dcd557`（22:46:57）。**收敛循环与对照没有运行**，`scripts/run_a17b.py` 从未执行。时间为 2026-10-04 Asia/Singapore，取自 `date` 与 git。

- 用户决定：22:32（回答后紧接的 `date` 为 22:32:52）把 `BUDGET.json` 的单次求解上限提到 600 s、总墙钟提到 5400 s；22:57（紧接的 `date` 为 22:57:48）在 G1 未通过后选择暂停 SC-01。
- 边界前的检查（22:15–22:29，结果已写进边界前言）：成本探测，等价命令 `python scripts/a17b_g1_diagnostics.py cost M D ORDER`；5 GHz 网格减半诊断，`python scripts/a17b_g1_diagnostics.py rate`；对 (1,4) 二阶的一次 cProfile 剖析（临时脚本，未入库）。当时用的是暂存区里的临时脚本，`scripts/a17b_g1_diagnostics.py` 是之后按同样的计算整理入库的。
- G1（22:36–22:47）：
  1. `python -m unittest discover -s tests`：第一次 44 个测试中 1 个失败，即 50 MHz 二阶的追踪测试（neff 差 6.1×10⁻³，要求 ≤10⁻³）。之后的移位敏感性诊断见边界第 8 节，等价命令 `python scripts/a17b_g1_diagnostics.py shift`。记录这一失败后，追踪测试改到 2 GHz（那里本征问题条件良好），22:46:47 重跑，44 个全部通过。
  2. `python scripts/check_evidence.py --replay --images ../tmp/pdfs/a03`（22:45:24–22:45:56）：通过，21 个可重放实验逐字节一致，A17 按设计只做静态核对。
  3. 在临时副本里重跑 `run_a17.py`（22:43:03–22:45:24，`solve()` 已改）：状态 `NUMERICAL_FAIL_NOT_CONVERGED`，两轮的 36 个阈值判定和 18 个判据 2b 标志都与登记的 A17 相同。
- 复核（22:58–23:02）：`python scripts/a17b_g1_diagnostics.py cost 1 1 2` 与 `shift` 各跑一次。50 MHz 二阶的 |Z0| 移位漂移为 0.92%（(4,1)）、1.9%（(4,1) 与 (2,1)），一阶为 1.3×10⁻⁴–5.0×10⁻⁴；500 MHz 二阶为 5.2×10⁻⁵（第一次 1.8×10⁻⁵）。数值每次不同，量级与第 8 节一致。
- 输入：本目录 `BUDGET.json`。参考点文件只读了频率列（边界前言），没有参与任何计算。
- 输出：无。没有 results.json，没有图。
- 环境：同 A17（Python 3.14.6、numpy 2.5.3、scipy 1.18.1、femwell `be2c547`、scikit-fem 12.0.2、gmsh 4.15.2；Apple M5、16 GB）。`run_a17b.py` 另需 psutil 7.2.2（已在 `requirements-lock.txt`）。
- 偏离：(1) G1 第 1 条未通过，按边界第 3 节在运行前停止（第 8 节）。(2) 追踪测试在记录失败后改到 2 GHz，原 50 MHz 的结论保留在测试说明与边界第 8 节。(3) `run_a17b.py` 的看门狗没有在真实运行中检验过。
