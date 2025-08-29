#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, os, re
from collections import defaultdict

# 3-letter AA -> 1-letter
AA3 = {
 "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
 "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
 "THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O"
}

# 人源偏好（简化版，每个氨基酸挑一个常见密码子）
HUMAN_CODON = {
 "A":"GCC","R":"CGT","N":"AAC","D":"GAC","C":"TGC","Q":"CAG","E":"GAG","G":"GGC",
 "H":"CAC","I":"ATC","L":"CTG","K":"AAG","M":"ATG","F":"TTC","P":"CCC","S":"AGC",
 "T":"ACC","W":"TGG","Y":"TAC","V":"GTG","U":"TGC","O":"TTT"  # U/O很少见，占位
}

def read_pdb_to_sequences(pdb_path):
    """按链提取主链顺序的氨基酸序列（基于 CA 原子），返回 dict: chain->AA序列"""
    chains = defaultdict(list)
    seen = set()
    with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            resn = line[17:20].strip()
            chain = line[21:22]
            resi = line[22:26].strip()
            key = (chain, resi)
            if key in seen:
                continue
            seen.add(key)
            aa = AA3.get(resn.upper(), "X")
            chains[chain].append(aa)
    # 合成序列
    seqs = {ch: "".join(a for a in aas if a != "X") for ch, aas in chains.items()}
    return seqs

def back_translate(aa_seq):
    dna = ["ATG"]  # 加起始 ATG（如果首位不是M，仍然加ATG作为ORF起始）
    for i, aa in enumerate(aa_seq):
        if i == 0 and aa == "M":
            # 已经是起始 Met，保留一个 ATG 即可
            continue
        dna.append(HUMAN_CODON.get(aa, "GCT"))  # 不识别的给个 A
    dna.append("TAA")  # 简单给一个终止密码子
    return "".join(dna)

def write_fasta(path, name, seq):
    with open(path, "w", encoding="utf-8") as f:
        f.write(">"+name+"\n")
        # wrap 60 nt/aa per line
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60]+"\n")

def main():
    ap = argparse.ArgumentParser(description="Extract AA from PDB and reverse-translate to ORF")
    ap.add_argument("--pdb", required=True, help="input PDB path")
    ap.add_argument("--out_protein", required=True, help="output AA FASTA")
    ap.add_argument("--out_orf", required=True, help="output ORF FASTA")
    ap.add_argument("--name", default="CLDN18_2-CAR_scfv", help="FASTA entry name")
    args = ap.parse_args()

    seqs = read_pdb_to_sequences(args.pdb)
    if not seqs:
        raise SystemExit("No sequences parsed from PDB (check file).")

    # 1) 输出每条链并拼接一个 scFv 序列（按链名排序拼接）
    chains_sorted = sorted(seqs.keys())
    aa_concat = "".join(seqs[ch] for ch in chains_sorted)

    # 写蛋白FASTA（拼接版，链名版本可自己扩展）
    write_fasta(args.out_protein, args.name+"_AA", aa_concat)

    # 2) 反向翻译为 ORF（DNA），再保存
    dna = back_translate(aa_concat)
    write_fasta(args.out_orf, args.name+"_ORF", dna)

    print("[OK] AA FASTA:", args.out_protein)
    print("[OK] ORF FASTA:", args.out_orf)
    print("Chains parsed:", ",".join(chains_sorted), " lengths:", [len(seqs[c]) for c in chains_sorted])

if __name__ == "__main__":
    main()
