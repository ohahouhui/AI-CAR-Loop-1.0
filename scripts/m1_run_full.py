# scripts/m1_run_full.py
import os, sys, pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
proc_dir = os.path.join(BASE, 'dataprocessed')
res_dir  = os.path.join(BASE, 'results')
fig_dir  = os.path.join(BASE, 'resultsfigures')
tab_dir  = os.path.join(BASE, 'resultstables')
os.makedirs(res_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(tab_dir, exist_ok=True)

EXPR = os.path.join(proc_dir, 'M1_expr_log2.tsv')
CLI  = os.path.join(proc_dir, 'M1_clinical.tsv')

TARGET_GENE = 'CLDN18'   # CLDN18.2 属于 CLDN18 基因，bulk级别以 symbol 汇总

def load_data():
    if not (os.path.exists(EXPR) and os.path.exists(CLI)):
        print('[错误] 缺少预处理文件，请先运行: python scripts/m1_preprocess.py')
        sys.exit(1)
    X = pd.read_csv(EXPR, sep='\t', index_col=0)   # genes x samples
    meta = pd.read_csv(CLI, sep='\t')
    # y: 1=Tumor, 0=Normal
    y = (meta['Group'].astype(str).str.lower()=='tumor').astype(int).values
    # 对齐列顺序
    X = X[meta['SampleID'].tolist()]
    return X, y, meta

def compute_tsi(X, meta, gene):
    # “Tumor/Normal 均值”定义的简易 TSI：mean(Tumor) / (mean(Tumor)+mean(Normal))
    sids = meta['SampleID'].tolist()
    if gene not in X.index:
        return np.nan
    g = X.loc[gene, sids]
    tumor = g[meta['Group'].str.lower()=='tumor'].values
    normal = g[meta['Group'].str.lower()=='normal'].values
    if len(tumor)==0 or len(normal)==0:
        return np.nan
    tmean, nmean = tumor.mean(), normal.mean()
    denom = tmean + nmean if (tmean + nmean)!=0 else np.nan
    return float(tmean/denom) if pd.notnull(denom) else np.nan

def compute_auc_for_gene(X, y, gene):
    if gene not in X.index: 
        return np.nan
    scores = X.loc[gene].values  # 以表达量作为打分（简单基线）
    try:
        return float(roc_auc_score(y, scores))
    except Exception:
        return np.nan

def kfold_auc(X, y, gene, k=5):
    if gene not in X.index: return np.nan
    skf = StratifiedKFold(n_splits=min(k, sum(y==0), sum(y==1), 5), shuffle=True, random_state=42)
    vals = []
    for tr, te in skf.split(X.T, y):
        scores = X.loc[gene].values[te]
        try:
            vals.append(roc_auc_score(y[te], scores))
        except:
            pass
    return float(np.mean(vals)) if len(vals)>0 else np.nan

def main():
    X, y, meta = load_data()

    # 目标基因（可扩展：你可以在这里放入候选列表做排名）
    genes = list(set([TARGET_GENE]) & set(X.index))
    if len(genes)==0:
        print(f'[警告] 在表达矩阵中未找到 {TARGET_GENE}，请检查 expression.tsv 的基因symbol。')
        sys.exit(0)

    rows = []
    for g in genes:
        tsi = compute_tsi(X, meta, g)
        auc5 = kfold_auc(X, y, g, k=5)
        rows.append({'gene': g, 'TSI': tsi, 'AUC_5fold': auc5})

    df = pd.DataFrame(rows).sort_values(['TSI','AUC_5fold'], ascending=False)
    out_csv = os.path.join(res_dir, 'M1_antigen_ranking.csv')
    df.to_csv(out_csv, index=False)
    with open(os.path.join(tab_dir, 'M1_metrics.txt'), 'w', encoding='utf-8') as f:
        f.write(f'Mean AUC (5-fold): {df["AUC_5fold"].mean():.4f}\n')
        f.write(f'Target gene {TARGET_GENE} TSI: {df["TSI"].iloc[0]:.4f}\n')

    print('[OK] M1 完成：')
    print(' -', out_csv)
    print(' -', os.path.join(tab_dir, 'M1_metrics.txt'))
    print('提示：后续可把候选基因列表扩展为 [CLDN18, GPC3, HER2, CD19, ...] 做排名。')

if __name__ == '__main__':
    main()
