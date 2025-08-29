# scripts/m1_preprocess.py
import os, sys, pandas as pd, numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
d_tcga = os.path.join(BASE, 'data', 'TCGA_STAD')
out_dir = os.path.join(BASE, 'dataprocessed')
os.makedirs(out_dir, exist_ok=True)

EXP = os.path.join(d_tcga, 'expression.tsv')
CLI = os.path.join(d_tcga, 'clinical.tsv')

def main():
    if not (os.path.exists(EXP) and os.path.exists(CLI)):
        print('[错误] 缺少 expression.tsv 或 clinical.tsv，请先按 m1_fetch_tcga.py 的提示放好文件。')
        sys.exit(1)

    expr = pd.read_csv(EXP, sep='\t')
    # 约定第一列为基因symbol
    gene_col = expr.columns[0]
    expr = expr.set_index(gene_col)

    cli = pd.read_csv(CLI, sep='\t')
    # 约定包含两列：SampleID, Group(值为Tumor/Normal)
    assert {'SampleID','Group'}.issubset(set(cli.columns)), 'clinical.tsv 必须包含 SampleID / Group 列'

    # 交集样本
    samples = [s for s in cli['SampleID'].tolist() if s in expr.columns]
    expr = expr[samples]
    cli = cli[cli['SampleID'].isin(samples)].reset_index(drop=True)

    # 简单log2 转换（避免0）
    expr = np.log2(expr + 1)

    # 输出标准化文件
    expr_out = os.path.join(out_dir, 'M1_expr_log2.tsv')
    cli_out  = os.path.join(out_dir, 'M1_clinical.tsv')
    expr.to_csv(expr_out, sep='\t')
    cli.to_csv(cli_out, sep='\t', index=False)

    print('[OK] 预处理完成：')
    print(' -', expr_out)
    print(' -', cli_out)

if __name__ == '__main__':
    main()
