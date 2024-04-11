# coding:utf-8 
'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-04-06 13:58:33
LastEditors: Ne0tea
LastEditTime: 2024-04-11 16:12:06
'''
import sys
import os
import logging
import statistics
import pandas as pd
import numpy as np
from calculate_synteny_length_among_phy import calculate_synteny_region_length

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
def out_phy_ratio_calculate(row,phy_dic):
    """
    计算row.name与Out phylogeny共线区域长度的函数。
    参数:
    - row: pandas.DataFrame的一行，代表某个实体的数据。
    - phy_dic: 字典，键为实体名称，值为该实体对应的物理标识。
    返回值:
    - out_ratio_sum: 计算得到的外部物理比例之和。
    """
    # if row.name=='Chr12_18XHB-57_centro':
    #     print(row)
    out_ratio=[]
    out_ratio_dic={}
    same_ratio=[]
    cur_phy=phy_dic[row.name]
    for i in row.index:
        if i == row.name:
            continue
        if cur_phy != phy_dic[i]:
            out_ratio.append((i,row[i]))
            if phy_dic[i] in out_ratio_dic:
                out_ratio_dic[phy_dic[i]].append((i,row[i]))
            else:
                out_ratio_dic[phy_dic[i]]=[(i,row[i])]
        else:
            if np.isnan(row[i]):
                same_ratio.append((i,0))
            else:
                same_ratio.append((i,row[i]))
    
    # return out_ratio_sum/len(out_ratio)
    max_out_tuple = max(out_ratio, key=lambda x: float(x[1]))
    max_same_tuple = max(same_ratio, key=lambda x: float(x[1]))
    # mean_same_tuple = sum(float(x[1]) for x in same_ratio)/len(same_ratio)
    mean_out_tuple=[]
    # print(out_ratio_dic)
    for x in out_ratio_dic:
        cur_list=out_ratio_dic[x]
        cur_median=statistics.median([float(x[1]) for x in cur_list])
        mean_out_tuple.append((x,cur_median))
    mean_out_tuple = max(mean_out_tuple,key=lambda x: float(x[1]))
    mean_same_tuple = statistics.median([float(x[1]) for x in same_ratio])
    return str(round(max_out_tuple[1],4))+'|'+str(max_out_tuple[0].split('_')[1])+'|'\
        +phy_dic[max_out_tuple[0]]+'|'+str(round(max_same_tuple[1],4))+'|'+str(round(mean_same_tuple,4))+\
            '|'+str(round(mean_out_tuple[1],4))+'|'+str(mean_out_tuple[0])

def get_chr_synteny_length_ratio(SG_iden_bed,centro_file,phy_file,Out_CR,min_identity=95):
    out_value={}
    Out_CR_df=pd.read_csv(Out_CR,header=0)
    for index,row in Out_CR_df.iterrows():
        out_value[row['Chr_sp']+'_centro']=round(float(row['isOut']),5)
    centro_dic={}
    with open(centro_file,'r') as cf:
        for i in cf:
            line=i.strip().split()
            centro_dic[line[0]+'_centro']=int(line[2])-int(line[1])
    '''
    read stained glass bed file
    sg_dic* store bed file line as key
    identity_dic* store high identity region as chain
    '''
    try:
        last_chr1,last_chr2='',''
        with open(SG_iden_bed, 'r') as mf:
            sg_dic={}
            identity_dic={}
            synteny_len_dic={}
            for i in mf:
                line = i.strip().split()
                chr=line[0].split('_')[0]
                if i.startswith('#') or line[0] == line[3]:
                    continue
                if (last_chr1!=line[0] or last_chr2!=line[3]) and sg_dic and identity_dic:
                    '''if versus different then output last versus'''
                    condition,last_chr1,synteny_len1,last_chr2,synteny_len2=calculate_synteny_region_length((last_chr1,last_chr2),identity_dic)
                    if condition:
                        if last_chr1 in synteny_len_dic:
                            if last_chr2 in synteny_len_dic[last_chr1]:
                                synteny_len_dic[last_chr1][last_chr2]=synteny_len_dic[last_chr1][last_chr2]+synteny_len1/centro_dic[last_chr1]
                            else:
                                synteny_len_dic[last_chr1][last_chr2]=synteny_len1/centro_dic[last_chr1]
                        else:
                            synteny_len_dic[last_chr1]={last_chr2:synteny_len1/centro_dic[last_chr1]}
                    sg_dic={}
                    identity_dic={}
                '''restore synteny pair over min_identity'''
                if float(line[6]) > min_identity:
                    if line[2] in identity_dic:
                        identity_dic[line[2]].append(line[5])
                    else:
                        identity_dic[line[2]] = [line[5]]
                loc = '_'.join(line[1:3]) + '_' + '_'.join(line[4:6])
                sg_dic[loc] = i
                last_chr1,last_chr2=line[0],line[3]
            '''output the final versus'''
            condition,last_chr1,synteny_len1,last_chr2,synteny_len2=calculate_synteny_region_length((last_chr1,last_chr2),identity_dic)
            if last_chr1 in synteny_len_dic:
                if last_chr2 in synteny_len_dic[last_chr1]:
                    synteny_len_dic[last_chr1][last_chr2]=synteny_len_dic[last_chr1][last_chr2]+synteny_len1/centro_dic[last_chr1]
                else:
                    synteny_len_dic[last_chr1][last_chr2]=synteny_len1/centro_dic[last_chr1]
            else:
                synteny_len_dic[last_chr1]={last_chr2:synteny_len1/centro_dic[last_chr1]}
    except IOError as e:
        logging.error(f"File error: {e}")
        os.sys.exit(1)

    phy_classfly={}
    phy_order=[]
    with open(phy_file,'r') as phyf:
        for i in phyf:
            line=i.strip().split()
            phy_classfly[chr+'_'+line[0]+'_centro']=line[2]
            phy_order.append(chr+'_'+line[0]+'_centro')
    '''tranvert dict to dataframe
    '{'Chr01_13-65_centro': {'Chr01_18WR-118_centro': 0.05069965524234435,
                             'Chr01_18XHB-57_centro': 0.05069965524234435}}'
            Chr01_13-65_centro  B  C   D
        Chr01_18WR-118_centro  1  4  7  10
        Chr01_18XHB-57_centro  2  5  8  11
        Chr01_9311_centro  3  6  9  12
    '''
    cur_series_df=pd.DataFrame(synteny_len_dic,index=phy_order)
    '''fill na synteny species with 0'''
    cur_series_df = cur_series_df.fillna(0)
    cur_series_df = cur_series_df.T

    cur_series_df = cur_series_df.reindex(index=phy_order, columns=phy_order)
    cur_series_df = cur_series_df.fillna(0)
    Ge_Xi=[x for x in phy_order if phy_classfly[x] in ['GJ','XI']]
    # print(cur_series_df.index)
    # cur_series_df = cur_series_df.loc[Ge_Xi,Ge_Xi]
    out_phy_value=cur_series_df.apply(lambda row: out_phy_ratio_calculate(row,phy_classfly), axis=1)
    '''reindex the output with sp'''
    out_phy_value = out_phy_value.loc[Ge_Xi]
    out_phy_value.name=chr
    chr_out_ratio_df=pd.DataFrame(out_phy_value)
    chr_out_ratio_df['isOut']=chr_out_ratio_df.index.map(lambda x: out_value[x])
    chr_out_ratio_df['CR_len']=chr_out_ratio_df.index.map(lambda x: centro_dic[x])
    # out_phy_value.index = out_phy_value.index.map(lambda x: x.split('_')[1])
    return chr_out_ratio_df

if __name__ == "__main__":
    SG_iden_bed=sys.argv[1]#Absulate path shuold be sorted by column 1 and column 4
    # SG_iden_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_pw_SG\Chr12_MH63_YZ-2.test.bed'
    centro_file=sys.argv[2]
    # centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    phy_file=sys.argv[3]
    # phy_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_phylogenic.txt'
    Out_CR=sys.argv[4]
    # Out_CR=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\Out_CRL.csv'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        filename='log_all_70_synteny_record',
                        filemode='a')
    ratio_df=get_chr_synteny_length_ratio(SG_iden_bed,centro_file,phy_file,Out_CR)
    # print(ratio_df)
    ratio_df.to_csv(SG_iden_bed+'.csv')