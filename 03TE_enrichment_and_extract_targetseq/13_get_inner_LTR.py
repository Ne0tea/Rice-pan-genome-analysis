'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-10 23:42:56
LastEditors: Ne0tea
LastEditTime: 2024-01-14 23:21:59
'''
import sys
import pandas as pd
import re
import os
import logging
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def is_empty_file(file_path):
    if os.path.exists(file_path) and os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
        return False
    else:
        return True

def main(file,fastafile):
    raw_data=pd.read_table(file,header=None,usecols=[0,1,2,3,4,5,10],names=['Chr','Start','End','unit_name','nosense','ord','label'])
    raw_data['length']=raw_data['End']-raw_data['Start']
    pattern=re.compile(r'\"Motif:(.*)\"')
    grouped = raw_data.groupby('unit_name')
    pair_LTR=pd.DataFrame(columns=['Chr','Start','End','unit_name','nosense','ord','label'])
    # pair_LTR=pd.DataFrame()
    # 遍历每个分组
    inner_length=[]
    inner_unit_name=[]
    for group_name, group_df in grouped:
        # 找到当前分组中所有 'SZ-22_I' 元素的索引位置
        sz_22_i_indices = group_df.index[group_df['label'].str.contains('SZ-22_I')].tolist()
        above_sz_22_index = max(min(sz_22_i_indices) - 1, min(group_df.index))
        below_sz_22_index = min(max(sz_22_i_indices) + 1, max(group_df.index))
        # print(group_df.loc[above_sz_22_index:min(sz_22_i_indices)-1])
        max_col2_index = group_df.loc[group_df.loc[above_sz_22_index:min(sz_22_i_indices)-1]['length'].idxmax()].name
        max_col2_index2 =group_df.loc[group_df.loc[max(sz_22_i_indices)+1:below_sz_22_index]['length'].idxmax()].name
        nearest_indices=[max_col2_index2, max_col2_index]
        concensus_len=0
        count=0
        break_s=False
        c_LTR=pd.DataFrame(columns=['Chr','Start','End','unit_name','nosense','ord','label'])
        for i in nearest_indices:
            count+=1
            cur_len=group_df.loc[i]['End']-group_df.loc[i]['Start']
            if concensus_len !=0 and  abs(cur_len-concensus_len) > 0.1*concensus_len:
                break_s=True
                break
            else:
                concensus_len=cur_len
        if not break_s:
            # print((group_df.loc[sz_22_i_indices]['End']-group_df.loc[sz_22_i_indices]['Start']).values)
            inner_length.append((group_df.loc[sz_22_i_indices]['End']-group_df.loc[sz_22_i_indices]['Start']).values[0])
            inner_unit_name.append(group_df.loc[i]['unit_name'])
    num=0
    outfile=open(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_pairLTR\SZ-22_inner.txt','w')
    for i in range(0,len(inner_length)):
        line=inner_unit_name[i]+'\t'+'SZ-22_pairLTR'+str(num)+'\t'+str(inner_length[i])+'\t'+str(1-inner_length[i]/2912)
        outfile.write(line+'\n')
        # print(line)
        num+=1
    outfile.close()

if __name__ == "__main__":
    # annofile=sys.argv[1]#SZ-22_unit_centro.bed 
    annofile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_unit_centro.bed'
    # fastafile=sys.argv[2]#
    fastafile=r''
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(annofile,fastafile)