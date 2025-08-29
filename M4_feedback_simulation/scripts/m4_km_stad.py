#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

def read_table_any(path):
    return pd.read_csv(path, sep="\t", header=0, index_col=0, compression="infer")

def tcga_barcode15(x): 
    return str(x)[:15]

def pick_surv_cols(df):
    cols = {c.lower(): c for c in df.columns}
    os_flag = cols.get("os") or cols.get("overall_survival") or cols.get("event")
    os_time = cols.get("os.time") or cols.get("os_days") or cols.get("days_to_death") or cols.get("days_to_last_followup") or cols.get("days_to_last_follow_up")
    vital   = cols.get("vital_status")
    return os_flag, os_time, vital

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--gene", default="ENSG00000066405")  # CLDN18
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 读取表达矩阵
    expr = read_table_any(args.expr)
    expr.index = expr.index.to_series().astype(str).str.replace(r"\.\d+$", "", regex=True)
    if args.gene not in expr.index:
        raise SystemExit(f"{args.gene} not in expression matrix")
    g = expr.loc[args.gene].copy()
    g.index = [tcga_barcode15(c) for c in g.index]
    g = np.log1p(g)
    g = g.groupby(level=0).mean()

    # 读取表型数据
    ph = read_table_any(args.pheno).copy()
    if "sample" in ph.columns: 
        ph["barcode15"] = ph["sample"].map(tcga_barcode15)
    elif "submitter_id" in ph.columns: 
        ph["barcode15"] = ph["submitter_id"].map(tcga_barcode15)
    else:
        ph = ph.reset_index().rename(columns={"index": "sample"})
        ph["barcode15"] = ph["sample"].map(tcga_barcode15)

    # 生存列提取
    os_flag, os_time, vital = pick_surv_cols(ph)
    if os_flag and os_time in ph.columns:
        ph["event"] = pd.to_numeric(ph[os_flag], errors="coerce")
        ph["time"]  = pd.to_numeric(ph[os_time], errors="coerce")
    else:
        vs = ph.get(vital, pd.Series(index=ph.index, dtype=object)).astype(str).str.upper()
        d1 = pd.to_numeric(ph.get("days_to_death", np.nan), errors="coerce")
        d2 = pd.to_numeric(ph.get("days_to_last_followup", ph.get("days_to_last_follow_up", np.nan)), errors="coerce")
        ph["time"] = d1.fillna(d2)
        ph["event"] = (vs == "DEAD").astype(float)

    # 决定要合并的列
    need_cols = ["barcode15", "event", "time"]
    if "sample_type" in ph.columns:
        need_cols.append("sample_type")

    # 合并
    df = pd.merge(
        pd.DataFrame({"expr": g}),
        ph[need_cols],
        left_index=True,
        right_on="barcode15",
        how="inner"
    )

    # 缺失值处理 + Primary Tumor 过滤（如有）
    df = df.dropna(subset=["event", "time"])
    if "sample_type" in df.columns:
        df = df[df["sample_type"].str.contains("Primary Tumor", na=False)]

    # 分组
    cut = df["expr"].median()
    df["group"] = np.where(df["expr"] >= cut, "CLDN18-high", "CLDN18-low")

    # KM 曲线
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(5.5, 4))
    for label, color in [("CLDN18-high", "tab:red"), ("CLDN18-low", "tab:blue")]:
        sub = df[df["group"] == label]
        kmf.fit(sub["time"], sub["event"], label=label)
        kmf.plot(ci_show=False, color=color)
    a = df[df["group"] == "CLDN18-high"]
    b = df[df["group"] == "CLDN18-low"]
    p = logrank_test(a["time"], b["time"], event_observed_A=a["event"], event_observed_B=b["event"]).p_value
    plt.title(f"STAD OS by {args.gene} (median split)\nlog-rank p={p:.3g}")
    plt.xlabel("Days")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "M4_KM_STAD_CLDN18.png"), dpi=160)
    plt.close()

    # Cox 回归
    cph = CoxPHFitter()
    cdf = df[["time", "event", "expr"]].copy()
    cph.fit(cdf, duration_col="time", event_col="event")
    cph.summary.to_csv(os.path.join(args.outdir, "M4_Cox_summary.tsv"), sep="\t")

    # Meta 信息
    with open(os.path.join(args.outdir, "M4_KM_meta.txt"), "w") as f:
        f.write(f"gene={args.gene}\ncutoff=median={cut:.6g}\nN={len(df)} p_logrank={p:.6g}\n")

    print("[OK] KM + Cox done.")

if __name__ == "__main__":
    main()
