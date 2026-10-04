# 工具与日常工作方式

主线入口是 README → 立项书 → BOUNDARIES（v2）→ 来源覆盖矩阵 → 实验索引 → 当前实验的 BOUNDARY/INPUTS/RESULTS。

## 会使用什么工具

| 工具 | 用途 | 当前需要你做什么 |
|---|---|---|
| Python | 运行解析模型与验证 | 已装；不需要先系统学完整语言 |
| NumPy | 复数相位、频率数组和参数扫描 | 能看懂输入参数即可 |
| SciPy | 独立数值积分；后续插值/求根 | 已装，用来交叉验证解析式 |
| Matplotlib | 画平均电压、EO 频响、S11 与误差图 | 已装 |
| pandas | 保存和检查参数表、曲线数据 | 已装 |
| JupyterLab | 逐格运行学习笔记和小实验 | 已装，辅助理解 |
| unittest | 自动边界测试 | Python 自带；每个模块完成后运行 |
| Git / GitHub | 本地版本、公开展示、保存每一步证据 | 已有；当前任务会同步一次 |
| GitHub CLI（gh） | 创建仓库、推送状态检查 | 已登录，用于本次同步 |
| PDF 阅读工具 | 对照原文公式、坐标轴和图例 | 本地论文已有，通常直接打开即可 |
| WebPlotDigitizer 或同类工具 | 后续读取 Figure 2/3 曲线 | 到数字化阶段再决定，不需要现在安装 |
| VS Code | 可选代码编辑器 | JupyterLab 已足够入门 |
| COMSOL / HFSS / Lumerical | 几何和全波场求解 | 当前解析复现不需要 |

## 文件夹逻辑

```text
configs/                         参数及来源
src/tfln_mzm/                    可替换的物理模块
experiments/
  track_A_reproduction/          论文主线：每步边界、输出、结论
  track_B_extensions/            长度、损耗、速度失配扩展
tests/                           自动边界检查
data/digitized/                  论文曲线读取值与误差，待建立
data/simulated/                  辅助模拟数据
docs/                            边界、符号、来源、学习分配
  learning/                     课件静态背景与验证
notebooks/                       交互学习
scripts/                         一条命令重跑结果
references/                      引文与源文件指纹
```

代码按物理模块分，实验按研究问题分。不要为每个实验复制整套物理代码；通过独立配置调用同一核心，版本由 Git 记录。历史结果不因为新实验出现就被覆盖。

## 每一步怎么留证据

按 BOUNDARIES 第 7 节的 G0–G4 逐级放行（2026-10-04 起）：

1. G0：按 `experiments/EXPERIMENT_TEMPLATE.md` 写 BOUNDARY.md，回答运行前五问，判据带数字，**单独提交**后再运行。
2. G1：`python -m unittest discover -s tests` 与 `python scripts/check_evidence.py --replay` 都通过，并先复现父实验的关键数字。
3. G2/G3：运行模型，检查收敛、有限性、被动性和支持区间；写 INPUTS.md（命令、哈希、环境、偏离）。
4. G4：写 RESULTS.md，包括失败和不能证明的内容；把新输出登记进 `experiments/EVIDENCE_MANIFEST.json`。
5. 本地提交。推送到 GitHub 是对外发布，须先经用户同意（立项书第 0 节）。

论文 PDF、课件 PDF、虚拟环境和个人课程作业留在本地。公开仓库只放原创代码、说明、引用及允许分享的数据。后续数字化曲线同时保留来源和不确定度。
