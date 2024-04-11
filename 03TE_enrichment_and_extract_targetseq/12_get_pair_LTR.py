'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-10 23:42:56
LastEditors: Ne0tea
LastEditTime: 2024-02-27 23:40:16
'''
import sys
import pandas as pd
import re
import os
import logging
pd.options.display.max_rows = None
def is_empty_file(file_path):
    if os.path.exists(file_path) and os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
        return False
    else:
        return True

def main(file,fastafile,cluster_version:bool,TE='SZ-22'):
    if not cluster_version:
        raw_data=pd.read_table(file,header=None,usecols=[0,3,4,5,8,9,10],names=['Chr','unit_name','nosense','ord','Start','End','label'])
        raw_data.drop(raw_data[~(raw_data['label'].str.contains('Motif:SZ-22'))].index, inplace=True)
        raw_data.sort_values(by = ['Chr','Start'],ignore_index=True,inplace=True)
        #print(raw_data)
        raw_data['length']=raw_data['End']-raw_data['Start']
        grouped = raw_data.groupby('unit_name')
        pair_LTR=pd.DataFrame(columns=['Chr','Start','End','unit_name','nosense','ord','label'])
        # pair_LTR=pd.DataFrame()
        for group_name, group_df in grouped:
            sz_22_i_indices = group_df.index[group_df['label'].str.contains(TE+'_I')].tolist()
            max_col2_index = group_df.loc[group_df.loc[min(group_df.index):min(sz_22_i_indices)-1]['length'].idxmax()].name
            max_col2_index2 =group_df.loc[group_df.loc[max(sz_22_i_indices)+1:max(group_df.index)]['length'].idxmax()].name
            nearest_indices=[max_col2_index2, max_col2_index]
            concensus_len=0
            count=0
            break_s=False

            c_LTR=pd.DataFrame(columns=['Chr','Start','End','unit_name','nosense','ord','label'])
            for i in nearest_indices:
                count+=1
                cur_len=group_df.loc[i]['End']-group_df.loc[i]['Start']
                if concensus_len != 0 and  abs(cur_len-concensus_len) > 0.1*concensus_len:
                    break_s=True
                    break
                else:
                    concensus_len=cur_len
                c_line=pd.Series(group_df.loc[i]).to_frame().T
                c_line['unit_name']=c_line['unit_name']+'_'+str(count)
                c_LTR=pd.concat([c_LTR,c_line],ignore_index=True)
            if not break_s:
                pair_LTR=pd.concat([pair_LTR,c_LTR],ignore_index=True)

    if cluster_version:
        cmd_getbed='cut -f -4 '+file+' > tmp_pairLTR.bed'
        os.system(cmd_getbed)
    else:
        pair_LTR.to_csv('tmp_pairLTR.bed',sep='\t',index=0,header=0)

    logging.info('LTR anno achieve')
    LTR_dic={}
    os.system('mkdir -p pairLTR_bed')
    with open('tmp_pairLTR.bed','r') as bf:
        for i in bf:
            last_underscore_index = i.strip().split()[3].rfind('_')
            unit_name = i.strip().split()[3][:last_underscore_index]
            if unit_name in LTR_dic:
                LTR_dic[unit_name].append(i)
            else:
                LTR_dic[unit_name]=[i]
    for unit in LTR_dic:
        pairbedfile=open('./pairLTR_bed/'+unit+'.bed','w')
        for i in LTR_dic[unit]:
            pairbedfile.write(i)
        pairbedfile.close()
        cmd_getfa='bedtools getfasta -s -name -fi '+fastafile+' -bed '+'./pairLTR_bed/'+unit+'.bed'+' -fo '+unit+'.fa'
        os.system(cmd_getfa)
        os.system('sed -i \'s/::.*//g\' '+unit+'.fa')
        cmd_alifa='mafft '+unit+'.fa'+' > '+unit+'.fa_ali'
        os.system(cmd_alifa)
        os.system('rm '+unit+'.fa')
if __name__ == "__main__":
    # annofile=sys.argv[1]#SZ-22_unit_centro.bed 
    # annofile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_pairLTR\merge_full_LTR.bed'
    annofile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_pairLTR\SZ-22_unit_centro.bed'
    # fastafile=sys.argv[2]#
    fastafile=''
    # TE=sys.argv[3] # TE name ex, SZ-22
    # TE='SZ-22'

    # if sys.argv[3] == 0 :
    #     cluster_version=False
    # else:
    #     cluster_version=True
    cluster_version=False
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(annofile,fastafile,cluster_version)
