# M03 实际输入与运行记录

边界提交 `ab0497c`；M04 边界 `d432c27`（在本实验运行前）；脚本提交 `c8f1b39`（运行前）。运行于 2026-10-05 11:30（Asia/Singapore），墙钟 1.3 s。

- 命令：`python scripts/digitize_m03.py --images ../tmp/pdfs/a03`（橙、红、蓝阈值从 `scripts/digitize_a03.py` 导入；`to_phys`、`groups_of` 从 `scripts/digitize_m01.py` 导入）。
- 输入（只在本地，不进 Git）：`p4_1.png`（SHA-256 `6efa75979dca4e92`；与 `p4_0.jpeg` 同一次 `pdfimages -all -f 4 -l 4` 提取，见 M01 INPUTS）；PDF（`ddcf40b8514cd43d`）。`m03_metadata.json` 记录完整哈希。
- 输出：`data/digitized/liu2022/m03_curves.csv`（`41d8577b9fd0fa05`）、`m03_metadata.json`（`b1ff8ccb0e836231`）、本目录 `digitized_preview.png`（`ad937f3355888d77`）、`results.json`（`030c2adf658e75e8`）；已登记 `experiments/EVIDENCE_MANIFEST.json`。叠加图 `fig7_overlay.png` 写到本地栅格目录（含原图，不进 Git）：青色为判成虚线的像素，洋红/绿色为读数点。
- 确定性：连续运行两次，四个输出逐字节相同。
- G1：运行前 50 个测试通过；`check_evidence.py --replay --images ../tmp/pdfs/a03` 通过（静态 27 个登记实验，59 个输出逐字节相同）。运行前单独试过刻度与参考虚线的检测函数（只碰框线、刻度与灰黑虚线），没有运行曲线提取。
- 环境：Python 3.14.6、numpy 2.5.3、scipy 1.18.1、pandas 3.0.6、matplotlib 3.11.2、Pillow 12.3.0。
- 偏离：无。
