# extract_cldn18_counts.py
# 目的：用 Ensembl 基因ID 直接在 TCGA-STAD.star_counts.tsv 中提取 CLDN18 的表达
# 说明：Xena 的 STAR 矩阵通常是 log2(count+1) 值；第一列为 Ensembl_ID（含版本号）
import pandas as pd
from pathlib import Path

MATRIX = "TCGA-STAD.star_counts.tsv"
GENE_ID = "ENSG00000066405"  # CLDN18（基因层）

out_tables = Path("../resultstables")
out_tables.mkdir(parents=True, exist_ok=True)

print("[INFO] Loading matrix, this can take some seconds...")
df = pd.read_csv(MATRIX, sep="\t", index_col=0)  # 行：Ensembl_ID(可能含.版本)；列：样本

# 去掉行索引中的版本号（如 ENSGxxxx.xx -> ENSGxxxx）
idx_stripped = df.index.to_series().astype(str).str.split(".").str[0]
mask = idx_stripped.eq(GENE_ID)

if not mask.any():
    msg = f"CLDN18 ({GENE_ID}) not found in {MATRIX}"
    print("[WARN]", msg)
    (out_tables / "M1_CLDN18_gene_check.txt").write_text(msg, encoding="utf-8")
else:
    sub = df.loc[mask]
    # 若同一基因ID出现多行（不同版本），取行均值作为基因层表达
    gene_series = sub.mean(axis=0)
    desc = gene_series.describe()

    # 在全基因中计算变异度（方差/标准差分位）
    all_sd = df.std(axis=1)
    cldn18_sd = gene_series.std()
    sd_rank = int((all_sd > cldn18_sd).sum() + 1)
    sd_pct = 100.0 * sd_rank / len(all_sd)

    # 导出每样本的表达
    per_sample_csv = out_tables / "M1_CLDN18_counts_per_sample.csv"
    gene_series.to_csv(per_sample_csv, header=["log2_count_plus1"])

    # 导出文本报告
    report_lines = [
        f"Found rows (with version): {list(sub.index)}",
        f"Samples: {gene_series.shape[0]}",
        "Expression summary (log2(count+1)):",
        str(desc),
        f"Across-gene SD of CLDN18: {cldn18_sd:.4f}",
        f"SD rank among all genes: {sd_rank} / {len(all_sd)} (~{sd_pct:.2f} percentile)",
        f"Per-sample values saved to: {per_sample_csv.name}",
    ]
    report_txt = out_tables / "M1_CLDN18_gene_check.txt"
    report_txt.write_text("\n".join(report_lines), encoding="utf-8")

    print("[OK] Wrote", report_txt)
    print("\n".join(report_lines))
