# scripts/m1_fetch_tcia.py
import os
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
d = os.path.join(BASE, 'data', 'TCIA')

def main():
    os.makedirs(d, exist_ok=True)
    print('[TCIA] 可选文件（影像不参与AUC/TSI核心计算，仅用于M4联动）：')
    print(f' - {os.path.join(d,"tcia_manifest.csv")}  # 影像清单或索引')
    print('[提示] 若暂无，可先跳过，不影响M1全流程。')

if __name__ == '__main__':
    main()
