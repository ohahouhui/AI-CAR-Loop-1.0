# scripts/m1_fetch_spatial.py
import os
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
d = os.path.join(BASE, 'data', 'SPATIAL')

def main():
    os.makedirs(d, exist_ok=True)
    print('[Spatial] 可选文件（10x Visium等）：')
    print(f' - {os.path.join(d,"sp_counts.tsv")}   # spotx基因 矩阵')
    print(f' - {os.path.join(d,"sp_meta.tsv")}     # spot元数据（坐标/组织区域）')
    print('[提示] 若暂无，可先跳过。')

if __name__ == '__main__':
    main()
