 AI CAR Loop 1.0

A Modular AI-Driven Platform for Accelerated CAR-T Therapy Development, Illustrated with CLDN18.2-Positive Gastric Cancer

---Environment

OS: Windows 10/11
Python: 3.13.5 (venv recommended)
Dependencies:

  ```bash
  pip install pandas numpy matplotlib lifelines
  ```
Folder structure** (root = `C:\Users\surface\Desktop\AI-CAR-Loop-1.0`):

  ```
  AI-CAR-Loop-1.0/
    data/
    M1_antigen_discovery/
    M2_structure_docking/
    M3_mRNA_design/
    M4_feedback_simulation/
    M5_RL_loop/
    tank_out/
    scripts/
  

## M1 — High-Throughput Antigen Discovery (TANK)

**TANK** = **Target discovery by Analysis of geNe expression ranKing**
A variance-based high-throughput ranking method to identify candidate targets from RNA-seq.

### Data

* **Source**: UCSC Xena (TCGA Hub)
* **Cohort**: TCGA-STAD (Stomach Adenocarcinoma)
* **File**: `TCGA-STAD.star_counts.tsv.gz` (RNA-Seq STAR-counts, rows = genes, columns = samples)
* **Filtered genes**: 47,662 after min-detection filtering

### Command

```bash
cd C:\Users\surface\Desktop\AI-CAR-Loop-1.0\M1_antigen_discovery
python tank_rank.py ^
  --expr "C:\Users\surface\Desktop\AI-CAR-Loop-1.0\data\TCGA-STAD.star_counts.tsv.gz" ^
  --outdir "C:\Users\surface\Desktop\AI-CAR-Loop-1.0\tank_out" ^
  --targets ENSG00000066405 ENSG00000141736 ENSG00000120217
```

Where:

* `ENSG00000066405` = **CLDN18**
* `ENSG00000141736` = **ERBB2**
* `ENSG00000120217` = **CD274**

### Outputs

* `TANK_ranked.tsv` — full variance ranking
* `TANK_top100.tsv` — Top-100 ranked genes
* `TANK_targets.tsv` — rankings for specified targets
* `README_targets.txt` — run parameters + top-10 list

**Example run result**:

* CLDN18: Rank = 95 / 47,662
* ERBB2: Rank = 5,306
* CD274: Rank = 10,372

**Transition**: CLDN18 shows high expression heterogeneity; since CLDN18.2 is a clinically hot gastric cancer CAR-T target, we proceed to structural modeling in **M2**.

---

## M2 — Structural Modeling & Molecular Docking (CLDN18.2 case)

### Tools

* **AlphaFold2** Multimer mode (`model_2_multimer_v3`) — Google Colab
* **HADDOCK 2.4** Guru mode — web server

### Inputs

1. **scFv** (VH–linker–VL, from 14G11 mAb)
2. **Full CAR** (scFv + hinge + TM + costimulatory module 4-1BB/CD28 + CD3ζ)
3. **CLDN18.2** (UniProt Q8N6F1-2, full-length model from AlphaFold Protein Structure Database)

### AlphaFold2 Steps

* Build 3 inputs: scFv, full CAR, CLDN18.2
* Run each 5× independently
* Select model with pLDDT > 85 + optimal domain packing

### HADDOCK Docking Rounds

1. **Baseline** (scFv ↔ CLDN18.2 full) — BSA = 2524.1 Å², HADDOCK = -101.1, Z = -2.4
2. **Full CAR** (CAR full ↔ CLDN18.2 full) — Electrostatic = -135.3, HADDOCK = -108.8
3. **Refined re-dock** (scFv from round 2 ↔ CLDN18.2 fragment) — High energy, low Z (-0.7)
4. **ECL2 loop** (scFv ↔ CLDN18.2 ECL2) — Min RMSD, min restraint violations (11.7)

**Decision**: Docking #2 chosen for **M3**; Docking #4 kept as epitope-specific control.

---

## M3 — mRNA Design & Delivery Simulation

### Inputs

* `docking2_fullCAR_centroid.pdb` — from M2
* ORF from CLDN18.2–CAR (reverse-engineered from model)
* Platforms: **LNP**, **TMAB3**, **RNACap**

### Steps

1. **mRNA optimization** — codon optimize (human bias), GC \~ 55%, avoid >6bp repeats
2. **Delivery simulation** — Monte Carlo, 100 iterations, parameters:

   * Penetration
   * Selectivity
   * Stability
     Baselines: LNP 0.65±0.05, TMAB3 0.70±0.06 (+15% selectivity), RNACap 0.60±0.07 (GI-specific)

### Output example (top-5):

| Rank | Platform | Penetration | Selectivity | Stability | Score  |
| ---- | -------- | ----------- | ----------- | --------- | ------ |
| 1    | TMAB3    | 0.8632      | 1.15        | 0.88      | 0.8736 |
| 4    | RNACap   | 0.8697      | 1.05        | 0.92      | 0.8401 |

---

## M4 — In-Silico Feedback Simulation (CLDN18 Safety)

### Data

* Expression: `TCGA-STAD.star_counts.tsv.gz`
* Survival: `TCGA-STAD_curated_survival.txt`
* Gene: `ENSG00000066405` (CLDN18)
* Phenotype: TCGA barcode type (01 = Tumor, 11 = Normal)

### KM + Cox Command

```bash
cd C:\Users\surface\Desktop\AI-CAR-Loop-1.0\M4_feedback_simulation\scripts
python m4_km_stad.py ^
  --expr  "C:\Users\surface\Desktop\AI-CAR-Loop-1.0\data\TCGA-STAD.star_counts.tsv.gz" ^
  --pheno "C:\Users\surface\Desktop\AI-CAR-Loop-1.0\M4_feedback_simulation\input\TCGA-STAD_curated_survival.txt" ^
  --gene  "ENSG00000066405" ^
  --outdir "C:\Users\surface\Desktop\AI-CAR-Loop-1.0\M4_feedback_simulation\out"
```

### Safety Boxplot Command

```bash
python m4_safety_boxplot.py ^
  --expr  "C:\Users\surface\Desktop\AI-CAR-Loop-1.0\data\TCGA-STAD.star_counts.tsv.gz" ^
  --gene  "ENSG00000066405" ^
  --outdir "C:\Users\surface\Desktop\AI-CAR-Loop-1.0\M4_feedback_simulation\out"
```

Example results:

KM log-rank p = 0.882, HR = 1.341 (NS)
Tumor median > Normal median


## M5 — Reinforcement Learning Feedback Loop (Concept)

Purpose: show closed-loop from M4 → M1 using simulated RL.

Example table:

| Iter | min\_detect\_prop | Sim p\_logrank | Tumor/Normal |
| ---- | ----------------- | -------------- | ------------ |
| 1    | 0.10              | 0.882          | 1.85         |
| 3    | 0.20              | 0.520          | 2.35         |

---

Do you want me to also **embed the exact code snippets for each script** inside this README so PeerJ reviewers can run without opening extra files?
That would make it **100% self-contained** and even stronger for acceptance.
