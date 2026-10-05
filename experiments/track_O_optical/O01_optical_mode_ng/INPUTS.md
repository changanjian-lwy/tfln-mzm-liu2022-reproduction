# O01 实际输入与运行记录

边界提交 `65efa21`（12:40）；运行前修订 `62c741f`（13:02，BOUNDARY 第 5 节）；脚本提交 `3a3278e`（13:02，提交时没有求解过 Liu 截面）。运行于 2026-10-05 13:02:45–13:13:39（Asia/Singapore），墙钟 652.7 s，67 次求解，单次最长 47.4 s，峰值内存 1.40 GB；没有触发 `BUDGET.json` 的任何上限。

- G1 的非 FEM 部分（运行前）：62 个测试通过（50 个旧测试 + 12 个新测试）；`check_evidence.py` 静态通过（29 个登记实验）。
- 命令：`PYTHONPATH=src MPLBACKEND=Agg .venv/bin/python scripts/run_o01.py`（退出码 0）。
- 输入：`data/external/optical_materials/` 七个材料文件（SHA-256 前 16 位：CLN e `aae97bf9e9550a18`、CLN o `d254398af5abdff9`、MgO e `ccc7d932d0edde48`、MgO o `9315e2cc112ce89e`、Malitson `c587670f397241ec`、Ghosh-o `c63dabfdfe9921c8`、Olmon-ev `be778621e6491fc4`）；`BUDGET.json`（`b769d332f9c52c90`）。只在本地、脚本不读：引文 [10] PDF（`d62a71b62b2cf431`，波长与占空比的来源）。
- 输出：`results.json`（`d960cd8d13956dcf`）、`run_log.json`（`28e2696bd6233068`）、`o01_ng.png`（`f0d63d7caf90beb0`）；已登记 `experiments/EVIDENCE_MANIFEST.json`。只登记静态哈希：整轮 10.9 分钟，超过边界第 4 节的 10 分钟，不登记重放。
- 网格：接受 m = 2（半域；N 4144 单元，N 无电极 3002 单元，d = 1.5 时 5331 单元）；对照用的 m = 1 为 13755 单元。G1 平板 2904 单元。
- 环境：Python 3.14.6、numpy 2.5.3、scipy 1.18.1、scikit-fem 12.0.2、gmsh 4.15.2、shapely 2.1.2、matplotlib 3.11.2、femwell `be2c547`（`requirements-fem-lock.txt`，没有新装依赖）。
- 偏离：无。运行前的改动都在 BOUNDARY 第 5 节，判据数字未变。
- 运行前披露：(1) femwell 冒烟测试（圆截面棒，HE11）；(2) 非 Liu 截面（膜厚 1.0、刻蚀 0.5、脊宽 2.0、间隙 6 µm）干跑两次，全域 8 分钟、半域 69 s，判据阈值临时放宽，输出写在仓库外，不登记；(3) 第二次干跑的图中看到该非 Liu 截面的 ng 约 2.25（已写入 BOUNDARY 第 5 节）。
