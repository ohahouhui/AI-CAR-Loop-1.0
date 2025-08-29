# scripts/m1_fetch_scrna.py
import os
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
d = os.path.join(BASE, 'data', 'SC')

def main():
    os.makedirs(d, exist_ok=True)
    print('[scRNA] 可选文件（用于空间/单细胞佐证，不是M1硬依赖）：')
    print(f' - {os.path.join(d,"sc_counts.tsv")}   # 细胞x基因 计数矩阵')
    print(f' - {os.path.join(d,"sc_meta.tsv")}     # 细胞元数据（样本/细胞类型等）')
    print('[提示] 若暂无，可先跳过。')

if __name__ == '__main__':
    main()
