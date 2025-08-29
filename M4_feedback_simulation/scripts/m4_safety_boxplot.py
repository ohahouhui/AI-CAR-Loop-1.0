#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, argparse, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_table_any(path):
    return pd.read_csv(path, sep="\t", header=0, index_col=0, compression="infer")

def strip_version(x): 
    s = str(x)
    return s.split(".")[0] if "." in s and s.startswith("ENSG") else s

def parse_sample_type_from_tcga_barcode(barcode):
    """
    直接从表达矩阵的列名（TCGA 条形码）推断样本类型：
    第四段的前两位数字：01=Primary Tumor, 11=Solid Tissue Normal
    例: TCGA-AB-1234-01A-... -> '01' -> Tumor
        TCGA-AB-1234-11A-... -> '11' -> Normal
    """
    b = str(barcode)
    parts = b.split("-")
    if len(parts) >= 4:
        seg = parts[3]  # like '01A' or '11A'
        code = re.sub(r"[^0-9]", "", seg)[:2]
        if code == "01":
            return "Primary Tumor"
        if code == "11":
            return "Solid Tissue Normal"
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", required=True, help="TCGA-STAD star_counts tsv.gz")
    ap.add_argument("--gene", default="ENSG00000066405", help="CLDN18 Ensembl ID")
    ap.add_argument("--outdir", required=True)
    # --pheno 参数保留兼容，但不强制使用
    ap.add_argument("--pheno", default="", help="(optional) phenotype file; not required")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 读表达矩阵
    expr = read_table_any(args.expr)
    # 处理基因 ID
    expr.index = expr.index.to_series().map(strip_version)
    if args.gene not in expr.index:
        raise SystemExit(f"{args.gene} not found in expression matrix.")
    g = expr.loc[args.gene].copy()

    # 直接从列名解析样本类型
    types = [parse_sample_type_from_tcga_barcode(c) for c in g.index]
    df = pd.DataFrame({"sample": g.index, "expr": np.log1p(g.values), "sample_type": types})
    df = df[df["sample_type"].isin(["Primary Tumor", "Solid Tissue Normal"])]

    if df.empty:
        raise SystemExit("Could not infer any Primary Tumor / Solid Tissue Normal from expression column names. "
                         "请确认表达矩阵列名为标准 TCGA 条形码（如 TCGA-XX-XXXX-01A-...）。")

    # 出图
    plt.figure(figsize=(4.8, 4.2))
    df.boxplot(column="expr", by="sample_type", grid=False)
    plt.title(f"{args.gene} expression: Tumor vs Normal (log1p counts)")
    plt.suptitle("")
    plt.xlabel("")
    plt.ylabel("log1p(count)")
    outp = os.path.join(args.outdir, "M4_Safety_Tumor_vs_Normal.png")
    plt.tight_layout()
    plt.savefig(outp, dpi=160)
    plt.close()

    # 导出数值
    df[["sample","sample_type","expr"]].to_csv(os.path.join(args.outdir, "M4_Safety_values.tsv"), sep="\t", index=False)
    print("[OK] Safety plot saved:", outp, " N(Tumor)=", (df.sample_type=="Primary Tumor").sum(), 
          " N(Normal)=", (df.sample_type=="Solid Tissue Normal").sum())

if __name__ == "__main__":
    main()
