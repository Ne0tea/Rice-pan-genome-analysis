'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 19:49:56
LastEditors: Ne0tea
LastEditTime: 2024-01-24 21:49:28
'''
import sys
import re
import pandas as pd
import os

def main(anno_file,outfile):
    raw_data=pd.read_table(anno_file,header=None,names=['Chr','Unit_Start','Unit_End','unit_name','nosense','ord','Loc','Chr2','type_Start','type_End','label'])
    # raw_data['length']=raw_data['End']-raw_data['Start']
    pattern=re.compile(r'\"Motif:(.*)\"')
    grouped = raw_data.groupby('unit_name')
    result_df = pd.DataFrame(columns=raw_data.columns)
    count=0
    # 遍历每个分组
    for name, group in grouped:
        first_idx = group.index[0]
        last_idx = group.index[-1]
        print_statue=True
        for id in range(first_idx+1,last_idx):
            if 'LTR' in group.loc[id]['label'] and 'LTR' not in group.loc[id-1]['label'] and 'LTR' not in group.loc[id+1]['label']:
                print_statue=False
                print(group)
                count+=1
                break
        if print_statue:
            result_df = pd.concat([result_df, group])
    result_df.to_csv(outfile,sep='\t',index=False,header=False)
    print(count)
if __name__ == "__main__":
    #unit_anno_file=sys.argv[1]
    unit_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_unit_centroV2.detailed.bed'
    #filtered_file=sys.argv[2]
    filtered_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_unit_centroV2.detailed.filtered.bed'
    main(unit_anno_file,filtered_file)
