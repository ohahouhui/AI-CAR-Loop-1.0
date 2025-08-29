#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, sys, os, re, subprocess, shutil
from pathlib import Path

def read_fasta(p):
    if not p: return ""
    seq = []
    with open(p, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'): continue
            seq.append(line.strip().upper())
    return "".join(seq)

def gc_content(seq):
    if not seq: return 0.0
    g = seq.count('G'); c = seq.count('C')
    return (g+c)/len(seq)

def has_long_repeat(seq, k):
    if k <= 1: return False
    for base in "ACGUacguTt":
        b = base.upper().replace('T','U')
        if b*k in seq.replace('T','U'): return True
    return False

def check_mfe_with_tool(seq, tool="auto"):
    """尝试调用 RNAstructure 或 ViennaRNA（任选其一）；如果都没有则返回 None"""
    seq = seq.replace('T','U')
    tmpfa = "tmp_seq.fa"
    with open(tmpfa, "w") as f:
        f.write(">seq\n"+seq+"\n")
    mfe = None
    try:
        if tool in ("auto","vienna") and shutil.which("RNAfold"):
            # ViennaRNA
            out = subprocess.check_output(["RNAfold", "--noPS", tmpfa], text=True)
            for line in out.splitlines():
                if "(" in line and ")" in line:
                    # e.g.: .... (-23.40)
                    m = re.search(r"\(([-\d\.]+)\)", line)
                    if m: mfe = float(m.group(1))
        elif tool in ("auto","rnastructure") and shutil.which("Fold"):
            # RNAstructure CLI: Fold <seqfile> <ctfile>
            ct = "tmp.ct"
            subprocess.check_call(["Fold", tmpfa, ct], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            # 粗略估计：读取最后一行能量
            with open(ct, "r") as f:
                header = f.readline()
                m = re.search(r"ENERGY = ([-\d\.]+)", header.upper())
                if m: mfe = float(m.group(1))
            os.remove(ct)
    except Exception:
        pass
    finally:
        try: os.remove(tmpfa)
        except: pass
    return mfe

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orf", required=True)
    ap.add_argument("--utr5", default="")
    ap.add_argument("--utr3", default="")
    ap.add_argument("--target_gc", type=float, default=0.55)
    ap.add_argument("--avoid_repeats", type=int, default=6)
    ap.add_argument("--check_mfe", default="false")  # "true"/"false"
    ap.add_argument("--out", required=True)
    ap.add_argument("--lead_model", default="")  # 仅记录来源，非必需
    args = ap.parse_args()

    utr5 = read_fasta(args.utr5)
    orf  = read_fasta(args.orf)
    utr3 = read_fasta(args.utr3)
    if not orf:
        print("ERROR: ORF fasta is empty or not found.", file=sys.stderr); sys.exit(2)

    mrna = (utr5 + orf + utr3).replace("T","U")
    gc = gc_content(mrna)
    warn = []
    if gc < args.target_gc-0.1 or gc > args.target_gc+0.1:
        warn.append(f"GC {gc:.3f} deviates from target ~{args.target_gc:.2f}")
    if has_long_repeat(mrna, args.avoid_repeats):
        warn.append(f"Has >= {args.avoid_repeats} homopolymer run")

    mfe = None
    if str(args.check_mfe).lower() == "true":
        mfe = check_mfe_with_tool(mrna, tool="auto")
        if mfe is None:
            warn.append("MFE tool not found (RNAfold/Fold). Skipped ΔG check.")

    Path(os.path.dirname(args.out)).mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        f.write(">CLDN18_2-CAR_mRNA\n")
        f.write(mrna+"\n")

    # 同时写一份小日志
    meta = args.out.replace(".fasta", ".log.txt")
    with open(meta, "w") as f:
        f.write(f"lead_model={args.lead_model}\n")
        f.write(f"len={len(mrna)} gc={gc:.4f}\n")
        if mfe is not None:
            f.write(f"MFE_approx={mfe}\n")
        if warn:
            f.write("WARN="+" | ".join(warn)+"\n")
    print(f"[OK] Wrote {args.out}")
    if warn:
        print("[WARN]", " | ".join(warn))

if __name__ == "__main__":
    main()
