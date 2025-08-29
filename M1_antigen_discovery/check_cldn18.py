# check_cldn18.py
import pandas as pd
from pathlib import Path

MATRIX = "TCGA-STAD.star_counts.tsv"
PROBEMAP = "gencode.v36.annotation.gtf.gene.probemap"

out_dir = Path("resultstables")
out_dir.mkdir(parents=True, exist_ok=True)

# 读表达矩阵（行：EnsemblID.version，列：样本）
df = pd.read_csv(MATRIX, sep="\t", index_col=0)

# 读映射表（两列：id<tab>gene）
pm = pd.read_csv(PROBEMAP, sep="\t", header=None, names=["ensembl","gene"])
# 去掉版本号，增强匹配稳健性
pm["ensembl_stripped"] = pm["ensembl"].str.split(".").str[0]
df_index_stripped = df.index.to_series().str.split(".").str[0]

# 精确匹配 CLDN18（基因层）
ens_ids = pm.loc[pm["gene"]=="CLDN18","ensembl_stripped"].unique().tolist()
mask = df_index_stripped.isin(ens_ids)
sub = df.loc[mask]

report_lines = []
if sub.shape[0] == 0:
    report_lines.append("CLDN18 not found in this matrix at gene level.")
else:
    # 将多个转录变体合并为基因层（对多个行取列均值）
    gene_series = sub.mean(axis=0)
    desc = gene_series.describe()

    # 在全基因范围计算方差分位
    gene_sd = df.std(axis=1)
    cldn18_sd = gene_series.std()
    sd_rank = (gene_sd > cldn18_sd).sum() + 1
    sd_pct = 100.0 * (sd_rank / len(gene_sd))

    report_lines += [
        f"CLDN18 rows found: {list(sub.index)}",
        f"Samples: {gene_series.shape[0]}",
        "Expression (log2(count+1)) summary:",
        str(desc),
        f"Across-gene SD of CLDN18: {cldn18_sd:.4f}",
        f"SD rank among all genes: {sd_rank} / {len(gene_sd)} (~{sd_pct:.2f} percentile)",
    ]

out_path = out_dir / "M1_CLDN18_gene_check.txt"
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(report_lines))

print(f"[OK] Wrote {out_path}")
print("\n".join(report_lines))
