#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, json, os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def parse_pairlist(arg):
    # 形如: "LNP:0.65,0.05 TMAB3:0.70,0.06 RNACap:0.60,0.07"
    out = {}
    for token in arg.split():
        name, vals = token.split(":")
        mu, sd = vals.split(",")
        out[name] = (float(mu), float(sd))
    return out

def parse_scalars(arg):
    # 形如: "LNP:1.0 TMAB3:1.15 RNACap:1.05"
    out = {}
    for token in arg.split():
        name, val = token.split(":")
        out[name] = float(val)
    return out

def sample_truncnorm(mu, sd, n):
    x = np.random.normal(mu, sd, n)
    return np.clip(x, 0.0, 1.0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--platforms", default="LNP,TMAB3,RNACap")
    ap.add_argument("--priors", required=True)      # 均值,方差
    ap.add_argument("--selectivity", required=True) # 乘子
    ap.add_argument("--stability", required=True)   # 乘子
    ap.add_argument("--lead_model", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    np.random.seed(args.seed)
    platforms = [s.strip() for s in args.platforms.split(",") if s.strip()]
    priors = parse_pairlist(args.priors)
    sel    = parse_scalars(args.selectivity)
    stab   = parse_scalars(args.stability)

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    records = []
    for p in platforms:
        mu, sd = priors[p]
        pen = sample_truncnorm(mu, sd, args.iters)
        score = (pen * sel[p] * stab[p])  # 简单乘积，已截断0-1
        for i in range(args.iters):
            records.append({
                "platform": p,
                "penetration": float(pen[i]),
                "selectivity": float(sel[p]),
                "stability": float(stab[p]),
                "score": float(score[i])
            })

    # Top-5
    rec_sorted = sorted(records, key=lambda r: r["score"], reverse=True)
    top5 = rec_sorted[:5]
    with open(os.path.join(args.outdir, "delivery_top5.json"), "w") as f:
        json.dump({
            "lead_model": args.lead_model,
            "top5": top5
        }, f, indent=2)

    # 直方图
    plt.figure(figsize=(6,4))
    plt.hist([r["score"] for r in records], bins=20)
    plt.xlabel("Composite score"); plt.ylabel("Count"); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "delivery_hist.png"), dpi=160)
    plt.close()

    # 雷达图（平台均值对比）
    dims = ["penetration","selectivity","stability"]
    angles = np.linspace(0, 2*np.pi, len(dims), endpoint=False).tolist()
    angles += angles[:1]
    plt.figure(figsize=(5,5))
    ax = plt.subplot(111, polar=True)
    for p in platforms:
        mu, sd = priors[p]
        vals = [mu, sel[p], stab[p]]
        vals += vals[:1]
        ax.plot(angles, vals, label=p)
        ax.fill(angles, vals, alpha=0.1)
    ax.set_thetagrids(np.degrees(angles[:-1]), dims)
    ax.set_ylim(0, max(1.0, max(stab.values())*1.05))
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "delivery_radar.png"), dpi=160)
    plt.close()

    print("[OK] Wrote", os.path.join(args.outdir, "delivery_top5.json"))
    print("[OK] Figures:", "delivery_hist.png", "delivery_radar.png")

if __name__ == "__main__":
    main()
