#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TANK (Target discovery by Analysis of geNe expression ranKing)
Consolidated version (duplicate-index safe)
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np

# ===== Defaults (edit as needed) =====
DEFAULT_EXPR = r"C:\Users\surface\Desktop\AI-CAR-Loop-1.0\data\TCGA-STAD.star_counts.tsv.gz"
DEFAULT_OUTDIR = r"C:\Users\surface\Desktop\AI-CAR-Loop-1.0\tank_out"
DEFAULT_TARGETS = ["ENSG00000066405", "ENSG00000141736", "ENSG00000120217"]  # CLDN18, ERBB2, CD274
DEFAULT_IDTYPE = "auto"
DEFAULT_LOG1P = False
DEFAULT_MIN_DETECT_PROP = 0.10
DEFAULT_DETECT_THRESH = 1.0
DEFAULT_STAT = "var"          # var|mad|winsor
DEFAULT_WINSOR_ALPHA = 0.01
DEFAULT_TOPK = 100
DEFAULT_DUP_AGG = "none"      # mimic original behavior
# =====================================

def env_or_default(name, default, cast=str):
    val = os.environ.get(name, None)
    if val is None or val == "":
        return default
    try:
        return cast(val)
    except Exception:
        return default

def env_list(name, default_list):
    raw = os.environ.get(name, "")
    if not raw:
        return list(default_list)
    return [s.strip() for s in raw.split(",") if s.strip()] or list(default_list)

def infer_delimiter(path):
    with open(path, 'rb') as f:
        head = f.readline().decode('utf-8', errors='ignore')
    return '\t' if head.count('\t') >= head.count(',') else ','

def strip_ensembl_version_idx(idx):
    def _strip(x):
        if isinstance(x, str) and x.startswith('ENSG'):
            return x.split('.')[0]
        return x
    try:
        return idx.map(_strip)
    except Exception:
        return pd.Index([_strip(x) for x in idx])

def guess_id_type(index_like):
    sample = list(index_like)[: min(1000, len(index_like))]
    ens = sum(1 for g in sample if isinstance(g, str) and g.startswith('ENSG'))
    return 'ensembl' if ens > (len(sample) - ens) else 'symbol'

def winsorize_frame(df, alpha=0.01):
    alpha = float(alpha)
    if not (0.0 <= alpha < 0.5):
        raise ValueError("winsor_alpha must be in [0, 0.5)")
    qlow = df.quantile(alpha, axis=1)
    qhi  = df.quantile(1.0 - alpha, axis=1)
    return df.clip(lower=qlow, upper=qhi, axis=0)

def compute_score(df, stat='var', winsor_alpha=0.01):
    if stat == 'var':
        return df.var(axis=1, ddof=1)
    elif stat == 'mad':
        med = df.median(axis=1)
        return (df.sub(med, axis=0)).abs().median(axis=1)
    elif stat == 'winsor':
        wf = winsorize_frame(df, alpha=winsor_alpha)
        return wf.var(axis=1, ddof=1)
    else:
        raise ValueError("stat must be one of {'var','mad','winsor'}")

def load_listfile(path):
    if not path:
        return None
    items = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if s:
                items.append(s)
    return items

def load_targets(args, defaults):
    tlist = []
    if args.targets:
        tlist.extend(args.targets)
    if args.targets_file:
        tlist.extend(load_listfile(args.targets_file) or [])
    if not tlist:
        tlist = env_list("TANK_TARGETS", defaults)
    norm = [(t.split('.')[0] if isinstance(t, str) and t.startswith('ENSG') else t) for t in tlist]
    seen = set(); out = []
    for t in norm:
        if t not in seen:
            out.append(t); seen.add(t)
    return out

def drop_duplicates(df, how="none"):
    how = (how or "none").lower()
    if how == "none":
        return df
    if not df.index.duplicated(keep=False).any():
        return df
    if how == "mean":
        return df.groupby(df.index).mean()
    if how == "sum":
        return df.groupby(df.index).sum()
    if how == "first":
        return df[~df.index.duplicated(keep='first')]
    if how == "max":
        v = df.var(axis=1, ddof=1)
        tmp = v.reset_index()
        tmp.columns = ['gene', 'var']
        tmp['rowid'] = np.arange(len(tmp), dtype=int)
        pick = tmp.loc[tmp.groupby('gene')['var'].idxmax(), 'rowid'].values
        mask = np.zeros(len(df), dtype=bool)
        mask[pick] = True
        return df[mask]
    raise ValueError("dup_agg must be one of {'none','mean','sum','max','first'}")

def main():
    ap = argparse.ArgumentParser(description="TANK consolidated: variance-based ranking + target report (dup-safe)")
    ap.add_argument('--expr', help='Expression matrix path (.tsv/.csv/.gz)')
    ap.add_argument('--targets', nargs='+', help='Targets (same ID namespace as matrix)')
    ap.add_argument('--targets_file', help='File with targets (one per line)')
    ap.add_argument('--id_type', choices=['auto','ensembl','symbol'], help='ID namespace in expression index')
    ap.add_argument('--log1p', action='store_true', help='Apply log1p before scoring')
    ap.add_argument('--min_detect_prop', type=float, help='Keep genes detected in ≥ this fraction of samples (0 disables)')
    ap.add_argument('--detect_thresh', type=float, help='Detection threshold on ORIGINAL scale (0 disables)')
    ap.add_argument('--stat', choices=['var','mad','winsor'], help='Score statistic')
    ap.add_argument('--winsor_alpha', type=float, help='Winsor alpha (for --stat winsor)')
    ap.add_argument('--dup_agg', choices=['none','mean','sum','max','first'], help='Resolve duplicate gene IDs')
    ap.add_argument('--gene_list', help='Keep only genes listed in this file')
    ap.add_argument('--sample_keep', help='Keep only samples (columns) listed here')
    ap.add_argument('--topk', type=int, help='If >0, also write Top-K table')
    ap.add_argument('--outdir', help='Output directory')
    args = ap.parse_args()

    # Defaults + env
    expr = args.expr or env_or_default("TANK_EXPR", DEFAULT_EXPR, str)
    outdir = args.outdir or env_or_default("TANK_OUTDIR", DEFAULT_OUTDIR, str)
    id_type = args.id_type or env_or_default("TANK_IDTYPE", DEFAULT_IDTYPE, str)
    log1p = args.log1p or bool(int(env_or_default("TANK_LOG1P", "1" if DEFAULT_LOG1P else "0", int)))
    min_detect_prop = args.min_detect_prop if args.min_detect_prop is not None else float(env_or_default("TANK_MIN_DETECT_PROP", str(DEFAULT_MIN_DETECT_PROP), float))
    detect_thresh = args.detect_thresh if args.detect_thresh is not None else float(env_or_default("TANK_DETECT_THRESH", str(DEFAULT_DETECT_THRESH), float))
    stat = args.stat or env_or_default("TANK_STAT", DEFAULT_STAT, str)
    winsor_alpha = args.winsor_alpha if args.winsor_alpha is not None else float(env_or_default("TANK_WINSOR_ALPHA", str(DEFAULT_WINSOR_ALPHA), float))
    topk = args.topk if args.topk is not None else int(env_or_default("TANK_TOPK", str(DEFAULT_TOPK), int))
    dup_agg = args.dup_agg or env_or_default("TANK_DUP_AGG", DEFAULT_DUP_AGG, str)
    gene_list = args.gene_list or env_or_default("TANK_GENE_LIST", "", str)
    sample_keep = args.sample_keep or env_or_default("TANK_SAMPLE_KEEP", "", str)
    targets = load_targets(args, DEFAULT_TARGETS)

    if not expr or not os.path.exists(expr):
        sys.stderr.write("ERROR: Expression file not found.\n"
                         f"Provided/Default: {expr}\n")
        sys.exit(2)

    # Load matrix
    delim = infer_delimiter(expr)
    df = pd.read_csv(expr, sep=delim, header=0, index_col=0, compression='infer')
    df.index = strip_ensembl_version_idx(df.index)

    # Optional sample subset
    if sample_keep:
        keep_cols = set(load_listfile(sample_keep))
        exist = [c for c in df.columns if c in keep_cols]
        if not exist:
            sys.stderr.write("ERROR: No overlap between sample_keep and columns.\n"); sys.exit(3)
        df = df[exist]

    # Optional gene whitelist
    if gene_list:
        keep_rows = set(load_listfile(gene_list))
        keep_rows = {g.split('.')[0] if isinstance(g, str) and g.startswith('ENSG') else g for g in keep_rows}
        df = df.loc[df.index.intersection(keep_rows)]

    # ID namespace (informational)
    id_ns = guess_id_type(df.index) if id_type == 'auto' else id_type

    # Duplicate aggregation
    df = drop_duplicates(df, how=dup_agg)

    # Detection filter (compute base detect_prop if thresholds active)
    if detect_thresh > 0 or min_detect_prop > 0:
        detect_prop_base = (df > detect_thresh).sum(axis=1) / df.shape[1]
        keep = detect_prop_base >= float(min_detect_prop)
        df_filt = df.loc[keep].copy()
    else:
        df_filt = df

    # Transform
    if log1p:
        df_filt = np.log1p(df_filt)

    # Score + mean
    score = compute_score(df_filt, stat=stat, winsor_alpha=winsor_alpha)
    mean = df_filt.mean(axis=1)

    # IMPORTANT: compute detect_prop on the filtered matrix to avoid duplicate-index reindex
    if detect_thresh > 0 or min_detect_prop > 0:
        detect_prop_used = (df_filt > detect_thresh).sum(axis=1) / df_filt.shape[1]
    else:
        detect_prop_used = pd.Series(1.0, index=df_filt.index)

    ranked = pd.DataFrame({'score': score, 'mean': mean, 'detect_prop': detect_prop_used}).sort_values('score', ascending=False)

    # Outputs
    os.makedirs(outdir, exist_ok=True)
    ranked_path = os.path.join(outdir, 'TANK_ranked.tsv')
    ranked.to_csv(ranked_path, sep='\t')

    topk_path = None
    if topk and topk > 0:
        topk_path = os.path.join(outdir, f'TANK_top{topk}.tsv')
        ranked.head(topk).to_csv(topk_path, sep='\t')

    # Targets
    present_rows = []
    not_found = []
    norm_targets = [(t.split('.')[0] if isinstance(t, str) and t.startswith('ENSG') else t) for t in targets]
    for t in norm_targets:
        if t in ranked.index:
            rpos = int(ranked.index.get_loc(t)) + 1
            present_rows.append({
                'target': t,
                'rank': rpos,
                'score': ranked.loc[t, 'score'],
                'mean': ranked.loc[t, 'mean'],
                'detect_prop': ranked.loc[t, 'detect_prop']
            })
        else:
            not_found.append(t)
    targets_df = pd.DataFrame(present_rows).sort_values('rank') if present_rows else \
                 pd.DataFrame(columns=['target','rank','score','mean','detect_prop'])
    targets_path = os.path.join(outdir, 'TANK_targets.tsv')
    targets_df.to_csv(targets_path, sep='\t', index=False)

    # Report
    report_path = os.path.join(outdir, 'README_targets.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("TANK target ranking report (Consolidated, dup-safe)\n")
        f.write(f"Expression file: {expr}\n")
        f.write(f"Genes after dup_agg/gene_list: {df.shape[0]}\n")
        f.write(f"Genes after filter: {ranked.shape[0]}\n")
        f.write(f"Samples: {df.shape[1]}\n")
        f.write(f"ID namespace: {id_ns}\n")
        f.write(f"dup_agg: {dup_agg}\n")
        f.write(f"Preprocessing: log1p={'on' if log1p else 'off'}, min_detect_prop={min_detect_prop}, detect_thresh={detect_thresh}\n")
        f.write(f"Statistic: {stat} (winsor_alpha={winsor_alpha if stat=='winsor' else 'NA'})\n")
        if gene_list:
            f.write(f"Gene whitelist applied: {gene_list}\n")
        if sample_keep:
            f.write(f"Sample subset applied: {sample_keep}\n")
        f.write("\nTop-10 genes by score:\n")
        top10 = ranked.head(10).reset_index()
        for i, row in top10.iterrows():
            g = row.iloc[0]
            f.write(f"{i+1}. {g}\t score={row['score']:.6g}\t mean={row['mean']:.6g}\t detect_prop={row['detect_prop']:.3f}\n")
        f.write("\nRequested targets found:\n")
        if present_rows:
            for row in present_rows:
                f.write(f"  {row['target']}\t rank={row['rank']}\t score={row['score']:.6g}\t mean={row['mean']:.6g}\t detect_prop={row['detect_prop']:.3f}\n")
        else:
            f.write("  (none)\n")
        if not_found:
            f.write("\nRequested targets NOT FOUND in index (check ID namespace and spelling):\n")
            for t in not_found:
                f.write(f"  {t}\n")

    print("[TANK] Wrote:")
    print(" ", ranked_path)
    if topk_path:
        print(" ", topk_path)
    print(" ", targets_path)
    print(" ", report_path)

if __name__ == '__main__':
    main()
