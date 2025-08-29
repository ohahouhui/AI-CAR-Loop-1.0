import pandas as pd

# 读取文件（如果是gz压缩格式，可以直接读取）
file_path = "TCGA-STAD.star_counts.tsv.gz"

# 读取前几行看结构
df = pd.read_csv(file_path, sep="\t")

# 总行数（基因数）
total_genes = df.shape[0]

# 总列数（包含第一列基因ID）
total_cols = df.shape[1]

# 样本数（除去第一列）
total_samples = total_cols - 1

# 打印前5个样本ID
sample_ids = df.columns[1:6].tolist()

print("总基因数:", total_genes)
print("总样本数:", total_samples)
print("前5个样本ID:", sample_ids)
