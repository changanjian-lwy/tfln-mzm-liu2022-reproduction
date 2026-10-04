# 交接文件

更新：2026-10-04。写给下一个接手本项目的人（或 AI 会话）。先读本文，再读 [立项书](PROJECT_CHARTER.md) 和 [边界 v2](docs/BOUNDARIES.md)。

## 1. 现在在哪

- 复现对象：Liu et al., IEEE PTL 34(16), 854–857 (2022)，容性加载 TFLN 调制器，公式 (1)–(4) 与 Figure 3。
- 已完成：A00–A16（Track A）、B01–B04（Track B），共 21 个登记实验，全部可以逐字节重放。
- 一句话结论：只用论文 Fig.2 作输入、不拟合，40 Ω 的平均电压与 EO 响应、50 Ω 的 EO 响应（含首次 −3 dB 下穿，差 0.04 GHz）在读数精度内复现；S11 面板与开路负载没有复现。17 个可评分系列：4 个通过、5 个失败在输入读数精度内、8 个超出（其中 6 个未解释）。详见 [Figure 3 逐系列结论](docs/FIGURE3_STATUS.md)。
- 阶段：P0–P3 完成；P4（用户验收）进行中。SC-01 第一阶段已授权、未开始。

## 2. 怎么上手

```sh
python3 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements-lock.txt && python -m pip install -e . --no-deps
python -m unittest discover -s tests                 # 38 个测试
python scripts/check_evidence.py --replay            # 重放全部登记实验（A03/A06 需本地栅格，见第 6 节）
```

各实验入口：`scripts/run_a0*.py`、`scripts/run_a1*.py`、`scripts/run_b01.py`、`scripts/run_b02_b04.py --study B02|B03|B04`、`scripts/digitize_a03.py`/`digitize_a06.py`。

## 3. 规则（必须遵守）

1. **边界先行**：按 [实验模板](experiments/EXPERIMENT_TEMPLATE.md) 写 BOUNDARY.md（运行前五问、带数字的判据、预测与判据分开、停止条件），**单独提交**后再运行。
2. **放行关卡 G0–G4**：G1 要求测试与 `check_evidence.py --replay` 通过，并先复现父实验的关键数字。
3. **登记证据**：运行后写 INPUTS.md、RESULTS.md，用 `python scripts/register_evidence.py '<json>'` 登记。已登记文件只读；发现错误写进 [勘误](docs/ERRATA.md)，不改原文件。
4. **不拟合未知量**：不调 Z0、损耗、nm、f0 去对 Figure 3，不改图例，不填数字化缺口（边界第 12 节）。对不上就登记冲突分支，交给用户。
5. **模型层级**：每个实验只声明一层（L0 解析、L1 数字化、L2 FEM、L3 网络），只按该层评分。
6. **授权边界**：范围变更、放宽判据、推送、付费资源须先问用户（立项书第 0 节）。

## 4. 已做的决定与仍待决的事项

完整记录见 [决策记录](docs/DECISIONS_PENDING.md)。

| # | 事项 | 状态 |
|---|---|---|
| D1 | 合并并推送 | 已决定（2026-10-04）：合并到 `main` 并推送 |
| D2 | 主线 f0 | 已决定：保持 1 MHz；A11 的约 1 GHz 结果作为敏感性结论并列报告 |
| D3 | Fig.3(c) "40 Ω" 曲线 | 已决定：按图例原样评分，结论中并列说明 A10 |
| D5 | SC-01 截面 FEM | 已决定：只授权第一阶段（方法关卡）；第二阶段需先找到几何来源并向用户报告 |
| D4 | 与 Fig.2(b) 比较时的带宽定义 | 待决；建议主线保持首次下穿，比较时两种定义都报 |
| D6 | SC-02 馈线/焊盘网络 | 待决；建议在 SC-01 之后 |
| D7 | SC-03 数字化测量图 Fig.6/7 | 待决；可做，但独立评分 |

## 5. 下一步建议（按顺序）

1. **SC-01 第一阶段（A17）**：安装 femwell（GPL-3.0，作为依赖使用，不复制代码），按其 `RF_CPW_transmission_line_tutorial` 复现 Tuncer 1994 的微波折射率与损耗。BOUNDARY 先写：网格收敛、外边界距离收敛、与参考数据的事前容差；从这一阶段起 `BUDGET.json` 生效（边界第 10 节）。通过后才能把 femwell 当作可信工具。
2. **第二阶段的准备（不运行）**：论文给出了非加载段截面（金厚 4 µm、信号线 80 µm、间隙 20 µm、600 nm TFLN / 300 nm 刻蚀、100 nm SiO2、BCB 1.5 µm、石英衬底、2 µm 键合层），但**没有给 T 形电极尺寸与周期**。先查引文 [8]、[10]（同一课题组前作）能否提供；再列出假设清单向用户报告。容性加载是沿传播方向的周期结构，需要加载段与非加载段分别求解再级联，单个截面不够。
3. D4、D6、D7 由用户决定后再做。
4. 课程方面（不在本仓库）：EEK5103 Part II 已开课，CA2 于 2026-11-03 开始、11-09 截止。

## 6. 容易踩的坑

- **时间**：记录里的时间一律取 `date` 或 git 提交时间，不要估计（勘误 E-01 的教训）。登记前把 RESULTS 里引用的每个数字与 results.json 核对（E-02）。
- **本地栅格**：A03/A06 的重放需要论文原生栅格 `p2_2.jpeg`、`p3_0.png`（不进 Git），放在仓库上级目录 `../tmp/pdfs/a03/`，PDF 放在 `../tmp/pdfs/reference_tfln_lpt2022.pdf`。用 `--images ../tmp/pdfs/a03` 重放；没有栅格时这两项会被跳过，只做哈希检查。
- **脚本互相导入**：`run_a10` 导入 `run_a05`；`run_a11` 导入 `run_a08`、`run_a09`；`run_a16`、`run_b01` 导入 `run_a09`；`run_b02_b04` 导入 `run_a07`、`run_b01`。改其中任何一个，都要跑 `--replay` 确认依赖它的实验仍逐字节一致。
- **复数 Z0**：`termination.py` 从 A12 起显式支持 Re>0 的复数 Z0；以前会静默丢掉虚部。实数路径没有变。
- **Z0 读数缺口**：Fig.2(c) 在 Z0 接近 50 Ω 的峰附近有缺口（约 9.6–10.7、19.6–21.5、29.6–30.7、65–70、78–84、95–98.5、106–108、190.4–191.5、194.8–195.9 GHz）。指标只用受支持频点；缺口里的极值可能被漏掉（B01 已说明）。
- **带宽定义**：首次 −3 dB 下穿对纹波窄谷极其敏感（A07、B01）。报带宽时写明定义，截尾时写 `>fmax`。

## 7. 文件地图

| 位置 | 内容 |
|---|---|
| `PROJECT_CHARTER.md` | 立项书与授权记录 |
| `docs/BOUNDARIES.md` | 边界 v2：范围、冻结参数、层级、关卡、禁止事项 |
| `docs/SOURCE_COVERAGE_MATRIX.md` | 论文逐项状态与冲突分支 C-01–C-07 |
| `docs/ASSUMPTION_CROSSCHECK.md` | 论文与当前模型的假设对照 |
| `docs/FIGURE3_STATUS.md` | Figure 3 逐系列结论 |
| `docs/DECISIONS_PENDING.md` | 决策记录与待决事项 |
| `docs/WORKLOG.md`、`docs/ERRATA.md` | 工作日志、勘误 |
| `docs/learning/session_2026-10-04_zh.md` | 讲解提纲与自测题 |
| `experiments/README.md` | 实验索引 |
| `experiments/EVIDENCE_MANIFEST.json` | 登记证据的哈希清单 |
| `src/tfln_mzm/` | 物理模块（SI 单位） |
| `data/digitized/liu2022/` | A03/A06 读数与元数据 |
