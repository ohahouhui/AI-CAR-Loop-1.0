import pandas as pd
from pathlib import Path

# -------------------- 参数配置 --------------------
MATRIX = "TCGA-STAD.star_counts.tsv"
CLDN18_ID = "ENSG00000066405"  # CLDN18 基因 ID
TOP_N = 100  # 输出前 N 名标准差排名的基因

out_tables = Path("../resultstables")
out_tables.mkdir(parents=True, exist_ok=True)

# -------------------- 1. 读取表达矩阵 --------------------
print("[INFO] Loading matrix, please wait...")
df = pd.read_csv(MATRIX, sep="\t", index_col=0)

# 去除版本号（ENSG00000123456.1 -> ENSG00000123456）
stripped_idx = df.index.to_series().astype(str).str.split(".").str[0]
df.index = stripped_idx

# -------------------- 2. 计算每个基因的标准差 --------------------
gene_sd = df.std(axis=1).sort_values(ascending=False)

# -------------------- 3. 输出前 TOP_N 的基因 --------------------
top_sd_df = gene_sd.head(TOP_N).reset_index()
top_sd_df.columns = ["Gene_ID", "Standard_Deviation"]
top_file = out_tables / f"M1_Top{TOP_N}_SD_genes.csv"
top_sd_df.to_csv(top_file, index=False)

# -------------------- 4. 特别输出 CLDN18 的信息 --------------------
cldn18_expr = df.loc[CLDN18_ID] if CLDN18_ID in df.index else None

report_lines = []

if cldn18_expr is not None:
    cldn18_sd = cldn18_expr.std()
    cldn18_rank = int((gene_sd > cldn18_sd).sum()) + 1
    cldn18_pct = 100.0 * cldn18_rank / len(gene_sd)

    # 保存每个样本的表达
    cldn18_expr.to_csv(out_tables / "M1_CLDN18_counts_per_sample.csv", header=["log2_count_plus1"])

    # 汇总报告
    report_lines += [
        f"[INFO] CLDN18 ({CLDN18_ID}) found in matrix.",
        f"Standard deviation: {cldn18_sd:.4f}",
        f"Rank: {cldn18_rank} / {len(gene_sd)} (~{cldn18_pct:.2f} percentile)",
        "Expression summary:",
        str(cldn18_expr.describe())
    ]
else:
    report_lines += [f"[WARN] CLDN18 ({CLDN18_ID}) not found in expression matrix."]

# -------------------- 5. 保存报告 --------------------
report_file = out_tables / "M1_Gene_SD_Report.txt"
report_file.write_text("\n".join(report_lines), encoding="utf-8")
print("[OK] Wrote:", report_file)

# -------------------- 控制台输出 --------------------
print("\n".join(report_lines))
print(f"[OK] Top {TOP_N} gene SD written to: {top_file}")
