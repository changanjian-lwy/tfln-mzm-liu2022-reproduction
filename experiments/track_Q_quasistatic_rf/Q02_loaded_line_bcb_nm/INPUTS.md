# Q02 实际输入与运行记录

边界提交 `493c548`（16:42:52，在 Q01 运行前）；运行前修订与脚本 `ac451d9`（17:09，写在 Q01 登记之后，只改方法）；修正 `0561c98`（17:19，见偏离 1）；只改绘图 `ec43606`（见偏离 2）。登记的运行：2026-10-05 17:34:09–17:53:05（Asia/Singapore），脚本自己的预算计时 452.2 s（CPU 约 7.4 min；shell 计时 18 min 56 s，中间进程有一段没有运行，`perf_counter` 不计系统睡眠），266 次求解、178 个网格，单次求解最长 2.5 s，峰值内存约 1.1 GB，最大网格 146 482 单元；没有触发 `BUDGET.json` 的上限。

- G1 的非 FEM 部分（脚本提交前）：74 个测试通过；`check_evidence.py` 静态通过（32 个登记实验）。
- 命令：`PYTHONPATH=src .venv/bin/python scripts/run_q02.py`（退出码 0）。
- 输入：`data/digitized/liu2022/curves.csv`（A03，`ddf19d3509412656`，内电感上限同 Q01）；`data/digitized/liu2022/a06_curves.csv`（A06，`27c3c45e90d63dcc`，Fig.2(b) 读数）；Q01 `results.json`（`ea878e15ceee8717`，父实验 L_ext）；`data/external/rf_materials/SOURCE.md` 与 `dc_materials/SOURCE.md`（常数写在 `rf_quasistatic.py` 里，与两份一致）；`BUDGET.json`（`788bda55af64ea7e`）。只在本地、脚本不读：BCB 数据表 PDF（`fd75274f5d710e25`）。
- 输出：`results.json`（`e174ee8398d65c20`）、`q02_bcb_nm.png`（`ca544f011085aa43`）、`run_log.json`（只登记静态哈希）；已登记 `experiments/EVIDENCE_MANIFEST.json`，可重放（需要可选 FEM 依赖）。
- 网格：接受 m = 2。N 在 t = 1.5 µm：U 16 131、H 19 535、N 21 239 单元（m = 2）；H 为 66 852（m = 1）、43 669（d = 2）。
- 环境：同 Q01，没有新装依赖。

## 偏离
1. **第一次运行崩溃（17:10:11–17:19:12，退出码 1）。** 选角点去重时遇到"同某角点"的条目（KeyError），异常没有被捕获，没有写出也没有打印任何结果。修正 `0561c98` 只改这一行；方法与判据都没有变。
2. **第二次运行完整（17:19:30–17:26:22，`0561c98`），看过结果后只改了图。** 图 (b) 的图例挡住了一个角点；`ec43606` 只把图例移到右上角。另在临时副本里重跑一次（17:26:56–17:33:45）：`results.json` 与图逐字节相同。第三次运行作为登记运行，`results.json` 与第二次逐字节相同。
- 运行前披露：非 Liu 截面的干跑（见 BOUNDARY 第 5 节）。
