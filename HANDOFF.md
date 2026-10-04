# 交接文件

更新：2026-10-04 22:04。写给下一个接手本项目的人（或 AI 会话）。先读本文，再读 [立项书](PROJECT_CHARTER.md) 和 [边界 v2](docs/BOUNDARIES.md)。**下一步是 A17b，已获授权，可以直接开始（第 5 节）。**

## 1. 现在在哪

- 复现对象：Liu et al., IEEE PTL 34(16), 854–857 (2022)，容性加载 TFLN 调制器，公式 (1)–(4) 与 Figure 3。
- 已完成：A00–A17（Track A）、B01–B04（Track B），共 22 个登记实验。21 个可以逐字节重放；A17 只登记静态哈希，因为 FEM 求解每次运行有 10⁻⁴ 量级的差异（结论与全部判定可复现）。
- Figure 3 结论（L1，只用论文 Fig.2 作输入、不拟合）：40 Ω 的平均电压与 EO 响应、50 Ω 的 EO 响应（含首次 −3 dB 下穿，差 0.04 GHz）在读数精度内复现；S11 面板与开路负载没有复现。17 个可评分系列：4 个通过、5 个失败在输入读数精度内、8 个超出（其中 6 个未解释）。详见 [Figure 3 逐系列结论](docs/FIGURE3_STATUS.md)。
- 那 6 个未解释的失败，需要真实的 Z0(f)、nm(f) 才能继续，所以有 SC-01（截面 FEM）。它的第一步是方法关卡：先证明 femwell 能复现一个已发表的 CPW 基准（Tuncer 1994）。**A17 做了这一步，结果 `NUMERICAL_FAIL`**（见第 5 节），femwell 还不能当作 L2 工具。
- 阶段：P0–P3 完成；P4（用户验收）进行中。
- Git：A17 全部工作与本交接已于 2026-10-04 22:10 经用户同意推送到 `origin/main`。以后每次推送仍需用户同意。

## 2. 怎么上手

```sh
python3 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements-lock.txt && python -m pip install -e . --no-deps
# 可选的 FEM 层（A17 起需要；不装时 4 个 FEM 测试跳过）
python -m pip install -c requirements-lock.txt -r requirements-fem-lock.txt
python -m pip install --no-deps "femwell @ git+https://github.com/HelgeGehring/femwell.git@be2c54763f6215757fd2dac543c9385a82ca5b26"
python -m unittest discover -s tests                 # 42 个测试
python scripts/check_evidence.py --replay            # 重放登记实验（A03/A06 需本地栅格，见第 6 节；A17 跳过）
```

各实验入口：`scripts/run_a0*.py`、`scripts/run_a1*.py`、`scripts/run_b01.py`、`scripts/run_b02_b04.py --study B02|B03|B04`、`scripts/digitize_a03.py`/`digitize_a06.py`；A17 是 `scripts/run_a17.py`（约 3 分钟，按设计以退出码 2 结束）。

## 3. 规则（必须遵守）

1. **边界先行**：按 [实验模板](experiments/EXPERIMENT_TEMPLATE.md) 写 BOUNDARY.md（运行前五问、带数字的判据、预测与判据分开、停止条件），**单独提交**后再运行。运行前要改方法，就在 BOUNDARY 末尾追加带时间的修订并提交（A17 第 5 节是例子），判据数字不能动。
2. **放行关卡 G0–G4**：G1 要求测试与 `check_evidence.py --replay` 通过，并先复现父实验的关键数字。
3. **登记证据**：运行后写 INPUTS.md、RESULTS.md，用 `python scripts/register_evidence.py '<json>'` 登记。已登记文件只读；发现错误写进 [勘误](docs/ERRATA.md)，不改原文件。失败也要登记。
4. **不拟合未知量**：不调 Z0、损耗、nm、f0 去对 Figure 3，不改图例，不填数字化缺口（边界第 12 节）。对不上就登记冲突分支，交给用户。
5. **模型层级**：每个实验只声明一层（L0 解析、L1 数字化、L2 FEM、L3 网络），只按该层评分。
6. **授权边界**：范围变更、放宽判据、推送、付费资源须先问用户（立项书第 0 节）。FEM 实验必须有 `BUDGET.json`，由运行脚本强制执行。

## 4. 已做的决定与仍待决的事项

完整记录见 [决策记录](docs/DECISIONS_PENDING.md)。

| # | 事项 | 状态 |
|---|---|---|
| D1 | 合并并推送 | 已决定：合并到 `main` 并推送（`fde1a91`） |
| D2 | 主线 f0 | 已决定：保持 1 MHz；A11 的约 1 GHz 结果作为敏感性结论并列报告 |
| D3 | Fig.3(c) "40 Ω" 曲线 | 已决定：按图例原样评分，结论中并列说明 A10 |
| D5 | SC-01 截面 FEM | 已决定：只授权第一阶段（方法关卡）；第二阶段需先找到几何来源并向用户报告 |
| D8 | A17 失败后 SC-01 如何继续 | **已决定（22:04）**：开新尝试 A17b，判据不变，只改数值方法；A17b 仍不过就停下报告，建议暂停 SC-01 |
| D4 | 与 Fig.2(b) 比较时的带宽定义 | 待决；建议主线保持首次下穿，比较时两种定义都报 |
| D6 | SC-02 馈线/焊盘网络 | 待决；建议在 SC-01 之后 |
| D7 | SC-03 数字化测量图 Fig.6/7 | 待决；可做，但独立评分 |
| — | 推送 A17 与本交接 | 已决定（22:10）：用户同意，已推送 |

## 5. 下一步：A17b（已授权）

### 5.1 A17 发生了什么

见 [A17 结果](experiments/track_A_reproduction/A17_femwell_cpw_method_gate/RESULTS.md) 与 [运行记录](experiments/track_A_reproduction/A17_femwell_cpw_method_gate/INPUTS.md)。

- 两轮收敛循环都没通过，第三轮需要 m=0.125，超出事前下限 0.25，于是 `NUMERICAL_FAIL`。按规则没有进入对照，**也没有人看过模型与 Tuncer 参考点的比较**。
- 失败原因一：48 MHz 加最细网格（m=0.25、d=2，11.7 万单元）时本征问题崩溃，一致性检查差 13–15%，nm 跳 15%。
- 失败原因二：|Z0| 每次网格减半仍变 1.07–1.28%（判据 ≤1%），而且减小得很慢。
- 已经收敛的部分：1 GHz 与 38.9 GHz 的 nm、α 对网格和计算域都收敛；计算域 d=2 已够大。

### 5.2 A17b 怎么做

父实验 A17（`59d5811`）。层级 L2 方法关卡；等级上限 `DIGITIZED_COMPARISON_PASS`，范围仅限本基准。

**不能改的**（照抄 A17 边界）：
- 结构：信号 7 µm、间隙 10 µm、地 100 µm、厚 0.8 µm、σ=6×10⁷ S/m、εr=13。
- 参考点文件 `data/external/tuncer1994_via_femwell/reference_points.csv`。
- 收敛频率 0.0482、1.0、38.9 GHz；收敛循环规则：从 (1,1) 起，m 下限 0.25，d 上限 8。
- 判据 2（Re neff ≤0.5%、α ≤2%、|Z0| ≤1%）、2b、3、4 的全部阈值。
- `BUDGET.json` 的上限；i0 仍用传导电流。

**只改两处数值方法**，都要在运行前写进 BOUNDARY：
1. **低频模式追踪**：A17 用 femwell 默认移位（1.1·k0²·max Re εr）。A17b 改为：同一频率下，用已求解的较粗一级网格（没有时用较小一级计算域）的 neff 作 `compute_modes(n_guess=…)`；第一级 (1,1) 保持默认移位（A17 中它在三个频率都通过了 2b）。目的是防止 48 MHz 崩溃。崩溃的根因（选错模式，还是病态）**没有确认**。
2. **二阶单元**：`compute_modes(order=2)`，针对 |Z0| 收敛慢。建议选它而不是角点加密：它只改一个开关、不改几何；但角点场奇异可能仍然限制收敛速度，这是已知风险。最终选哪个，在 BOUNDARY 里写明理由。

**代码组织**：
- 给 `src/tfln_mzm/cpw_fem.py` 的 `solve()` 加 `order=1, n_guess=None` 两个关键字参数。默认值必须保持 A17 的行为。
- 新写 `scripts/run_a17b.py`，可以从 `run_a17` 导入 `Budget`、`checks` 等。**不要改 `run_a17.py`**，它是已登记证据的入口。

**G1 要先做到**：
- 新测试：order=2 时直流电阻测试仍在 2% 以内；复数 `n_guess` 能被 femwell 接受（未验证过）。
- 全部测试通过，`check_evidence.py --replay` 通过。
- 父实验关键数字：在临时副本里重跑 `run_a17.py`，状态与 36 个阈值判定要和登记的 A17 相同。

**预算**：A17 的一阶单元在 11–12 万单元时每次求解 20–25 s，进程峰值内存 3.0 GB（来自 `run_log.json`）。二阶单元的自由度约 3–4 倍，单次求解预计接近 180 s 的上限。这是估计，没有测过。触发上限就停下报告，不要自行提高上限。

**重放**：在边界里事先写明——A17b 只登记静态哈希，另在临时副本重跑一次，比较状态与全部阈值判定（同 A17）。

**停止条件**：
- A17b 通过 → 才能考虑第二阶段的准备（见 5.3），并向用户报告。
- A17b 不通过 → 不要再开 A17c，停下报告用户，建议暂停 SC-01。
- 任何情况下都不放宽阈值，不改 σ、εr、几何、参考点去凑结果。

### 5.3 之后（A17b 通过才做）

第二阶段的准备（不运行）：论文给出了非加载段截面（金厚 4 µm、信号线 80 µm、间隙 20 µm、600 nm TFLN / 300 nm 刻蚀、100 nm SiO2、BCB 1.5 µm、石英衬底、2 µm 键合层），但**没有给 T 形电极尺寸与周期**。先查引文 [8]、[10]（同一课题组前作）能否提供，再列出假设清单向用户报告。容性加载是沿传播方向的周期结构，需要加载段与非加载段分别求解再级联，单个截面不够。另外，A17 只覆盖到 39 GHz、0.8 µm 金属；Liu 需要到 200 GHz、4 µm 金和 LN 各向异性，每换一个截面都要重做收敛检查。

课程方面（不在本仓库）：EEK5103 Part II 已开课，CA2 于 2026-11-03 开始、11-09 截止。

## 6. 容易踩的坑

**通用**
- **时间**：记录里的时间一律取 `date`、git 提交时间或文件时间戳，不要估计（勘误 E-01）。登记前把 RESULTS 里引用的每个数字与 results.json 核对（E-02）。
- **本地栅格**：A03/A06 的重放需要论文原生栅格 `p2_2.jpeg`、`p3_0.png`（不进 Git），放在仓库上级目录 `../tmp/pdfs/a03/`，PDF 放在 `../tmp/pdfs/reference_tfln_lpt2022.pdf`。用 `--images ../tmp/pdfs/a03` 重放；没有栅格时这两项会被跳过，只做哈希检查。
- **脚本互相导入**：`run_a10` 导入 `run_a05`；`run_a11` 导入 `run_a08`、`run_a09`；`run_a16`、`run_b01` 导入 `run_a09`；`run_b02_b04` 导入 `run_a07`、`run_b01`。改其中任何一个，都要跑 `--replay` 确认依赖它的实验仍逐字节一致。
- **复数 Z0**：`termination.py` 从 A12 起显式支持 Re>0 的复数 Z0；以前会静默丢掉虚部。实数路径没有变。
- **Z0 读数缺口**：Fig.2(c) 在 Z0 接近 50 Ω 的峰附近有缺口（约 9.6–10.7、19.6–21.5、29.6–30.7、65–70、78–84、95–98.5、106–108、190.4–191.5、194.8–195.9 GHz）。指标只用受支持频点；缺口里的极值可能被漏掉（B01 已说明）。
- **带宽定义**：首次 −3 dB 下穿对纹波窄谷极其敏感（A07、B01）。报带宽时写明定义，截尾时写 `>fmax`。

**FEM（A17 起）**
- **安装**：femwell 必须 `--no-deps` 安装。它声明的依赖 meshwell 会拉进 cadquery、vtk、trame，求解用不到。装 FEM 依赖时加 `-c requirements-lock.txt`，确保原有版本不变。
- **电流定义**：教程的回路电流 ∮H·dl 不收敛（A17 边界第 5 节），模块用信号线内的传导电流 ∫σE_z dA。不要改回去。
- **低频下限**：本结构在 10 MHz 时本征问题给出无意义的结果；48 MHz 加最细网格时也崩溃过。每个求解点都要看判据 2b（γ 与 Z0 的 RLGC 一致性），不能只看 neff。
- **运行间噪声**：同一输入每次求解，neff 差到 1.4×10⁻⁴、|Z0| 差到 4.7×10⁻⁴（相对）。保留 6 位有效数字也做不到逐字节重放。
- **网格尺度**：m 只缩放网格尺寸，过渡距离不变；d 放大计算域，地宽取 min(100 µm, 域内)。
- **防止偷看**：基准结构上 0.05–39 GHz 的任何求解都接近某个参考点。G1 调试时只打印比值（R/R_dc、一致性），不打印 nm、α。真看到了，就在边界修订里披露（A17 就这样做过）。
- **退出码**：`run_a17.py` 状态不是 COMPLETE 时按设计以退出码 2 结束。
- **参考点来源**：28 个点由 femwell 作者数字化，是测量还是计算没有核实（见 `data/external/tuncer1994_via_femwell/SOURCE.md`）。

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
| `docs/environment.md`、`requirements-lock.txt`、`requirements-fem-lock.txt` | 环境与锁定版本（FEM 层单独锁定） |
| `docs/learning/session_2026-10-04_zh.md` | 讲解提纲与自测题 |
| `experiments/README.md` | 实验索引 |
| `experiments/EVIDENCE_MANIFEST.json` | 登记证据的哈希清单 |
| `src/tfln_mzm/` | 物理模块（SI 单位）；`cpw_fem.py` 是 FEM 截面求解（只在 L2 用） |
| `scripts/run_a17.py` | A17 运行脚本（收敛循环、预算强制、对照） |
| `data/digitized/liu2022/` | A03/A06 读数与元数据 |
| `data/external/tuncer1994_via_femwell/` | Tuncer 基准参考点与来源说明 |
