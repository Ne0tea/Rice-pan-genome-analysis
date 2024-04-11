'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-02-16 21:42:43
'''
import sys
import logging
import matplotlib.pyplot as plt
import numpy as np

def main(monomer_fai_file,LTR_insert_loc_file):
    monomer_len_dic={}
    ID_loc=[]
    with open(monomer_fai_file,'r') as mf:
        for i in mf:
            line=i.strip().split()
            monomer_len_dic[line[0]]=int(line[1])
    with open(LTR_insert_loc_file,'r') as LTRf:
        for i in LTRf:
            line=i.strip().split()
            # print(line)
            monomer_len=monomer_len_dic[line[1]]
            loc_ratio=round(float(line[2])/monomer_len,4)
            ID_loc.append(loc_ratio)
    print(min(ID_loc))
    # unique_values, frequencies = np.unique(ID_loc, return_counts=True)

    # 创建直方图
    plt.hist(ID_loc,bins=40, edgecolor='black')


    # 添加标题和标签
    plt.title('Frequency Histogram')
    plt.xlabel('Values')
    plt.ylabel('Frequency')

    # 显示图形
    plt.show()
if __name__ == "__main__":
    # monomer_fai_file=sys.argv[1]
    monomer_fai_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\Selected_70m_representative_protable.monomr.fasta.fai'
    # LTR_insert_loc_file=sys.argv[2]
    LTR_insert_loc_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_ID_monomer.txt'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(monomer_fai_file,LTR_insert_loc_file)
