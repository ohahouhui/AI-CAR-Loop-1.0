#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt
GENES = ["ENSG00000153563","ENSG00000172116","ENSG00000100479","ENSG00000180644"]  # CD8A/B,GZMB,PRF1
def read_table_any(p): return pd.read_csv(p, sep="\t", header=0, index_col=0, compression="infer")
def strip_ver(s): return str(s).split(".")[0]
def tcga_barcode15(x): return str(x)[:15]
ap = argparse.ArgumentParser()
ap.add_argument("--expr", required=True); ap.add_argument("--gene", default="ENSG00000066405"); ap.add_argument("--outdir", required=True)
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
expr = read_table_any(a.expr); expr.index = expr.index.to_series().map(strip_ver)
need = [a.gene] + GENES; 
for g in need:
    if g not in expr.index: raise SystemExit(f"Missing {g}")
sub = expr.loc[need].applymap(lambda x: np.log1p(x))
sub.columns = [tcga_barcode15(c) for c in sub.columns]; sub = sub.groupby(axis=1, level=0).mean()
cldn = sub.loc[a.gene]; imm = sub.loc[GENES]
z = imm.apply(lambda col: (col - imm.mean(axis=1))/imm.std(axis=1), axis=0)
proxy = z.mean(axis=0)
joined = pd.DataFrame({"CLDN18":cldn, "ImmuneProxy":proxy}).dropna()
rho = joined.corr(method="spearman").iloc[0,1]
plt.figure(figsize=(4.5,4)); plt.scatter(joined["CLDN18"], joined["ImmuneProxy"], s=12, alpha=0.6)
plt.xlabel("CLDN18 log1p"); plt.ylabel("CD8A/B+GZMB+PRF1 z-mean"); plt.title(f"Spearman rho={rho:.2f}")
plt.tight_layout(); plt.savefig(os.path.join(a.outdir,"M4_ImmuneProxy_scatter.png"), dpi=160); plt.close()
joined.to_csv(os.path.join(a.outdir,"M4_ImmuneProxy_values.tsv"), sep="\t"); print("[OK] Immune proxy done.")
