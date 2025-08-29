# scripts/m1_fetch_tcga.py
import os, sys

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
d = os.path.join(BASE, 'data', 'TCGA_STAD')

def main():
    print('[TCGA] 目标文件（请手动下载后放到此目录）：')
    print(f' - {os.path.join(d, "expression.tsv")}  # 基因x样本 的表达矩阵（建议TPM/FPKM，含基因symbol列）')
    print(f' - {os.path.join(d, "clinical.tsv")}    # 样本ID/分组(如Tumor/Normal)等临床信息')
    os.makedirs(d, exist_ok=True)
    miss = []
    for fname in ['expression.tsv','clinical.tsv']:
        fp = os.path.join(d, fname)
        if not os.path.exists(fp):
            miss.append(fp)
    if miss:
        print('\n[提示] 未检测到以下文件：')
        for m in miss: print('  -', m)
        print('\n[操作] 到 GDC/GEO 下载 STAD 表达与临床（或你已有的全量文件），重命名为上述文件名后放入该目录。')
        sys.exit(0)
    else:
        print('[TCGA] 文件已就绪。')

if __name__ == '__main__':
    main()
