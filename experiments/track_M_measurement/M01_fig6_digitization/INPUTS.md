# M01 实际输入与运行记录

边界提交 `3ef0605`；M02 边界 `2b80199` 在本实验运行前提交；脚本提交 `f364ea9`（运行前）。运行于 2026-10-05 10:21（Asia/Singapore），墙钟 1.7 s。

- 栅格准备（运行前）：`pdfimages -all -f 4 -l 4 reference_tfln_lpt2022.pdf p4`（poppler 26.09.0），得到第 4 页的原始 JPEG 流 `p4-000.jpg`，复制为本地栅格目录中的 `p4_0.jpeg`（与 A03 的 `p2_*`、`p3_*` 同一命名方式）；`p4-001.png`（Fig.7）同时复制为 `p4_1.png`，本实验不用。
- 命令：`python scripts/digitize_m01.py --images ../tmp/pdfs/a03`（颜色阈值从 `scripts/digitize_a03.py` 导入）。
- 输入（只在本地，不进 Git）：`p4_0.jpeg`（SHA-256 `79198a8026bc874f`）；PDF（`ddcf40b8514cd43d`，与 A03 相同）。`m01_metadata.json` 记录完整哈希。
- 输出：`data/digitized/liu2022/m01_curves.csv`（`cdbb25e45e5bc5e8`）、`m01_metadata.json`（`463f5d1891cf2f7b`）、本目录 `digitized_preview.png`（`4d2a1dbc0a8bf722`）、`results.json`（`aa9a63e4347e920f`）；已登记 `experiments/EVIDENCE_MANIFEST.json`。脚本另把叠加图 `fig6_overlay.png` 写到本地栅格目录（含原图，不进 Git）。
- 确定性：同一输入连续运行两次，四个输出逐字节相同。
- G1：运行前 50 个测试通过；`check_evidence.py --replay --images ../tmp/pdfs/a03` 通过（静态 25 个登记实验；52 个输出逐字节相同，A17/A17b 按设计跳过）。
- 环境：Python 3.14.6、numpy 2.5.3、pandas 3.0.6、matplotlib 3.11.2、Pillow 12.3.0。
- 偏离：无。`check_evidence.py` 只改了 `--images` 的帮助文字（提到 `p4_0.jpeg`），检查逻辑没有改。
