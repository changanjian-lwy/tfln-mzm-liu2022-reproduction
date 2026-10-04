# A17 实际输入与运行记录

边界提交 `83310ff`（21:23:07）；运行前修订 `eb4c692`（21:30:51）；模块、测试、脚本与 FEM 依赖锁定 `ec14442`（21:31:59）。运行于 2026-10-04 21:32:05–21:35:19（Asia/Singapore，取自 `date`）。

- 命令：`python scripts/run_a17.py`。退出码 2 表示状态不是 COMPLETE，本次状态为 `NUMERICAL_FAIL_NOT_CONVERGED`。G1：`python -m unittest discover -s tests`（42 个通过）与 `python scripts/check_evidence.py --replay --images ../tmp/pdfs/a03`（21 个已登记实验逐字节一致）。
- 输入：`experiments/track_A_reproduction/A17_femwell_cpw_method_gate/BUDGET.json`；`data/external/tuncer1994_via_femwell/reference_points.csv`（SHA-256 `4c044dbb53c14bea`）。这个文件只用来取收敛检查的最低、最高频率（0.0482、38.9 GHz）；参考值本身只在对照阶段使用，本次没有进入对照。
- 输出：本目录 `results.json`（SHA-256 `815cc34eb8d6af62`）与 `run_log.json`（计时与内存，SHA-256 `ba67a7d194ede920`）；已登记在 `experiments/EVIDENCE_MANIFEST.json`。没有画图：脚本只在进入对照时才画。
- 环境：Python 3.14.6、numpy 2.5.3、scipy 1.18.1；femwell 提交 `be2c547`（`--no-deps` 安装，包元数据版本号显示为 0.0.1）、scikit-fem 12.0.2、gmsh 4.15.2、pygmsh 7.1.17、meshio 5.3.5、shapely 2.1.2。见 `requirements-fem-lock.txt`。Apple M5、16 GB 内存。
- 资源：用时 193 s，18 次求解，单次最长 25.1 s，最大网格 117,031 单元，峰值内存 3.04 GB，都在 `BUDGET.json` 上限内。
- 偏离：(1) i0 改为传导电流，见边界第 5 节（运行前提交）。(2) 披露：G1 调试时看到 (1,1)、50 MHz 下 nm≈8.5，判据没有改动。(3) 重放方式与边界第 4 节的计划不同：results.json 无法逐字节重放（见下一条），A17 只登记静态哈希，`check_evidence.py --replay` 跳过它。
- 重放确定性（在临时副本中完整重跑一次，21:39:07 写出结果，用时 165 s，退出码同为 2）：状态、两轮的加密顺序、全部 36 个阈值判定（2 轮 × 2 种加密 × 3 个频率 × 3 个量）都相同，结论可复现。但原始数值每次运行都有差异：除崩溃那一次求解外，neff 最多差 1.4×10⁻⁴（相对），|Z0| 4.7×10⁻⁴，α 7×10⁻⁵，C 3.7×10⁻⁴。差异来源可能是 ARPACK 的随机起始向量，以及金属 σ/(ωε0) 达到 10¹⁰ 时移位求逆的舍入误差。收敛判据里的"变化量"是两个接近数之差，保留 6 位有效数字也吸收不了这种差异，所以字节不同。这一噪声底比判据 2 的阈值小一个数量级以上。
