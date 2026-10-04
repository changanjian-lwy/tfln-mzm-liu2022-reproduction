# A06 实际输入与运行记录

边界提交 `b84ed7a`；脚本与结果在其后的提交中加入。运行于 2026-10-04 约 15:05（Asia/Singapore）。

- 命令：`python scripts/digitize_a06.py --images ../tmp/pdfs/a03`（颜色阈值从 `scripts/digitize_a03.py` 导入，两者定义一致）。
- 输入（只在本地，不进 Git）：`p2_2.jpeg`（SHA-256 `1558d1694d0e525e`）、`p3_0.png`（`ba0d0b95b1ff9d31`），与 A03 同源；`a06_metadata.json` 记录完整哈希。
- 输出：`data/digitized/liu2022/a06_curves.csv`、`a06_metadata.json`、本目录 `digitized_preview.png`、`results.json`；哈希登记在 `experiments/EVIDENCE_MANIFEST.json`。脚本另把 `fig2b_overlay.png`、`fig3a_main_overlay.png` 写到本地栅格目录（含原图，不进 Git）。
- 确定性：同一输入连续运行两次，`a06_curves.csv` 逐字节相同。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2、Pillow 12.3.0。
- 偏离：第一次运行后在脚本末尾加入写 `results.json` 的代码（判据 2–4、6 的自动评估），再重跑；定标、掩码、阈值和生成 CSV 的代码都没有改。改动后连续两次运行的 `a06_curves.csv` 逐字节相同；首次运行的哈希没有单独记录。
