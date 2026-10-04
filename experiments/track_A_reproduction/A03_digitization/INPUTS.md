# A03 实际输入与复核记录

**2026-10-04 补写（P1 既有证据审查），不是当时的运行记录。** BOUNDARY 与 RESULTS 在提交 `8741206`（2026-09-18）同时进入 Git，事前声明无法由 Git 历史独立证明。

- 命令：`python scripts/digitize_a03.py --images ../tmp/pdfs/a03`。
- 输入（只在本地，不进 Git）：论文 PDF（SHA-256 `ddcf40b8514cd43d`）中提取的原生栅格 `p2_2.jpeg`（`1558d1694d0e525e`）与 `p3_0.png`（`ba0d0b95b1ff9d31`）；完整哈希记录在 `data/digitized/liu2022/metadata.json`。2026-10-04 复核时本地文件哈希与元数据一致。
- 代码：`scripts/digitize_a03.py` 与 `src/tfln_mzm/digitized.py` 最后修改于 `8741206`。
- 输出：`data/digitized/liu2022/curves.csv`（`ddf19d3509412656`）、`metadata.json`（`529c2e9c9b8984c9`）、`digitized_preview.png`（`36d9d369e8d87df8`）。
- 环境：见 `docs/environment.md`；需要 Pillow。
- 复核：2026-10-04 用本地栅格在临时副本中重跑，三个输出逐字节一致。没有本地栅格的人只能做静态哈希检查。

P1 分类：**可引用**，等级 `PARTIAL_EXTRACTION_COMPLETE`。需补证：Fig.3(a) 0–200 GHz 主图、Fig.2(b) 尚未数字化（见来源覆盖矩阵第 4 节）。
