# 光子集成电路的设计流程与本项目的位置

2026-10-05 写（决策 D11）。目的：弄清"电子 IC 用 Verilog，光子 IC 用什么"，并检查本项目的描述方法是否完备。第 4、5 节的事实来自实际安装和运行，第 1 节的商业工具信息来自厂商文档。

## 1. 一句话

光子 IC **没有** Verilog 那样"写行为、自动综合成电路"的语言。器件要一个个设计，再拼起来。描述分两层：

- **电路层**：每个器件是一个频域 S 参数块（紧凑模型），电路就是这些块的网表。作用相当于电子的 SPICE / Verilog-A。
- **物理层**：用 Python 写参数化单元（PCell），生成 GDSII 版图。

| 电子 IC | 光子 IC | 常见工具 |
|---|---|---|
| Verilog（数字 RTL，可综合） | 没有对应物 | — |
| SPICE / Verilog-A | S 参数紧凑模型 + 网表；光子 Verilog-A 模型可在 Spectre 里与驱动电路联合仿真 | Lumerical INTERCONNECT、CML Compiler；开源 SAX |
| 标准单元库 + PDK | 代工厂 PDK：PCell + 紧凑模型 + DRC 规则 | Synopsys OptoCompiler（Lumerical 已随 Ansys 并入 Synopsys，2025-07）、Luceda IPKISS；开源 gdsfactory |
| 版图 | GDSII，多为弯曲多边形，由 Python 生成 | gdsfactory、KLayout |
| TCAD / 寄生提取 | 器件级场求解：光模式、FDTD、RF 电极 | Lumerical MODE/FDTD、Tidy3D（付费云）、femwell（本项目 A17）、HFSS/CST |
| DRC/LVS | 有，但比电子 EDA 原始 | KLayout 规则脚本 |

## 2. 本项目在哪一层

Eq.(1)–(4) 是**一个器件（行波 MZM 的 RF/EO 响应）的紧凑模型**，处在 Verilog-A / S 参数那一层，只是用 Python 写的。TFLN 行波调制器的紧凑模型在业界也是这种传输线模型。差别在输入：业界的 Z0、nm、α 由几何尺寸经场求解算出；本项目从论文 Fig.2 读出（L1）。

## 3. 描述是否完备

| 层 | 状态 | 复现 Fig.3 | 设计优化 |
|---|---|---|---|
| 行波相互作用 + 终端（Eq.1–3） | 已实现、已验证（A01/A02；第 5 节另有 SAX 交叉核对） | 需要 | 需要 |
| RF 截面 Z0、nm、α(f) | 读数（L1）；FEM 暂停（SC-01） | 需要 | 关键：真实旋钮是几何，改几何会同时改 Z0/nm/α |
| T 形电极 → 等效均匀线 | 隐含假设；周期 ≪ 波长时成立 | 需要 | 需要 |
| 光模式 / Vπ·L | 只假设 Vπ∝1/L | 不需要（Fig.3 已归一化） | 需要 |
| 光损耗随长度 | 无 | 不需要 | 需要（B05 偏爱长线的原因之一） |
| 馈线、焊盘、电阻寄生 | 无（SC-02） | 可能是 S11 失败的原因 | 需要 |
| 工艺容差 | 无 | 不需要 | 需要（B06 的可行孤岛） |

结论：对复现目标，层级正确，缺的是输入和馈线网络（已登记）；对设计优化，至少缺光损耗、几何驱动的 RF 参数、工艺容差三层。在 L1 上做 B05/B06 式优化为时过早；BO 方法本身的验证（B05）不受影响。

## 4. 实际装上的工具（开源 TFLN PDK）

- 环境：课程目录下的独立环境，Python 3.12.15（uv 安装）。依赖用 `uv pip install --exclude-newer 2026-04-14` 解析，版本见 `requirements-pdk-lock.txt`：gdsfactory 9.40.1、SAX 0.16.14、JAX 0.9.2、KLayout 0.30.6。项目主环境（Python 3.14.6）不变。
- PDK：Luxtelligence `lxt_pdk_gf`，提交 `d16db6a`，MIT 许可。它的自带测试 120 项全部通过。用最新的 gdsfactory 9.51.0 时有 17 项 GDS 回归测试不过：XOR 显示只差 LN_SLAB 层上约 25 µm × 0.5 µm 的一条（12.5 µm²），T 形电极单元完全相同。参考 GDS 最后一次重生成是 2026-04-14，所以按这一天固定依赖。两个示例 notebook 运行无误。
- 读源码并运行后确认的事实：
  1. 调制器紧凑模型只有直流：`eo_phase_shifter`、`mzm_unbalanced` 输入波长与静态电压，没有 RF 频率轴。
  2. T 形电极单元 `trail_cpw` 与高速 MZM `mzm_unbalanced_high_speed` **没有任何电路模型**。本项目的 Eq.(1)–(4) 补的正是这一块。
  3. `trail_cpw` 默认几何（从 GDS 顶点核对）：信号电极 21 µm；T 形电极颈部 1.5 µm、头部 1.5 µm；两侧 T 头之间的间隙 4.0 µm（裸电极间隙 10 µm），2.5 µm 宽的 LN 脊居中；头长 44.7 µm、头间距 5.0 µm、颈宽 7.0 µm；周期 49.7 µm。由此估算周期结构截止频率 c/(2·2.25·49.7 µm) ≈ 1.3 THz，远高于 200 GHz，等效均匀线的描述在这类尺寸下成立。**这是另一家代工厂、另一个材料栈（lnoi400）的尺寸，不是 Liu 的**（Liu 印出的 T 形电极间隙是 3 µm）。将来恢复 SC-01 第二阶段，它只能作为 `external_reference` 假设分支。

## 5. 桥接演示（学习用，不是登记证据）

`scripts/pdk_sax_rf_bridge.py`（需要上述可选环境）把 A04 L1 的电极写成 SAX 网表：参考端口 → 均匀线（Z0(f)、γ(f)、6 mm）→ 负载，再与项目闭式 `termination.s11` 比较。在 6 891 个受支持频点上，20/40/50/80 Ω 四个负载的最大差 ≤4e-16。也就是说，本项目的 S11 与标准 S 参数网表法完全一致，模型可以直接接进这套流程。

同一脚本在线前插入一个 20 fF 的并联电容，演示 L3（SC-02）的网络怎样描述：只要在网表里多加 2 端口块即可。**20 fF 是随手取的示意值，没有来源**，结果（40 Ω 时 0–200 GHz 最大 S11 从 −16.75 dB 变到 −4.90 dB）只说明这类寄生值得认真建模，不能当证据。

## 6. 怎样复现

```bash
brew install uv && uv python install 3.12
uv venv --python 3.12 .venv-pdk && source .venv-pdk/bin/activate
git clone https://github.com/Luxtelligence/lxt_pdk_gf && cd lxt_pdk_gf && git checkout d16db6a
uv pip install --exclude-newer 2026-04-14 -e ".[dev]" && python install_tech.py && pytest -q   # 120 passed
cd <本仓库> && python scripts/pdk_sax_rf_bridge.py
```

`install_tech.py` 只是把 KLayout 工艺文件软链接到 `~/.klayout/tech/lnoi400`。
