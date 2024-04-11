'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-19 17:39:59
LastEditors: Ne0tea
LastEditTime: 2023-12-19 10:10:11
'''
import sys
import pandas as pd
import os
import scipy.stats as stats
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)

def main(seq_stat_file,out_centro_file):
    seq_df=pd.read_csv(seq_stat_file)
    seq_df.columns = ['chr_sp'] + seq_df.columns[1:].tolist()
    OC_seq_df=pd.read_csv(out_centro_file)
    OC_seq_df.columns = ['chr_sp'] + OC_seq_df.columns[1:].tolist()

    seq_df[['chr','sp']]=seq_df['chr_sp'].str.split('_',expand=True)
    seq_df=seq_df.drop('chr_sp',axis=1)
    OC_seq_df[['chr','sp']]=OC_seq_df['chr_sp'].str.split('_',expand=True)
    OC_seq_df=OC_seq_df.drop('chr_sp',axis=1)
    # sp_group_df=seq_df.drop('chr',axis=1)
    # os.makedirs('parse_seq_step_tmp',exist_ok=True)
    # print(seq_df)
    sp_group_df=seq_df.groupby(['sp']).apply(sum).drop(columns=['chr','sp'], axis=1)
    chr_group_df=seq_df.groupby(['chr']).apply(sum).drop(columns=['chr','sp'], axis=1)
    OC_sp_group_df=OC_seq_df.groupby(['sp']).apply(sum).drop(columns=['chr','sp'], axis=1)
    OC_chr_group_df=OC_seq_df.groupby(['chr']).apply(sum).drop(columns=['chr','sp'], axis=1)
    sp_group_unit_dom_out=open(seq_stat_file.replace('01','02').split('.')[0]+'_dom_sp_long.txt','w')
    sp_group_unit_ful_out=open(seq_stat_file.replace('01','02').split('.')[0]+'_ful_sp_long.txt','w')
    sp_group_unit_dom_out.write('Sp\tname\tstatue\tvalue\tunitnumber\tOCvalue\tOCunitvalue\tcondition'+'\n')
    sp_group_unit_ful_out.write('Sp\tname\tstatue\tvalue\tunitnumber\tOCvalue\tOCunitvalue\tcondition'+'\n')
    for index, row in sp_group_df.iterrows():
        sp=index
        sp_unit_num=int(sp_group_df.loc[sp,'unit_number'])
        sp_out_unit_num=int(OC_sp_group_df.loc[sp,'unit_number'])
        for column, value in row.items():
            if column=='unit_number':
                continue
            seq=column
            if seq in OC_sp_group_df:
                OC_value=int(OC_sp_group_df.loc[sp,seq])
            else:
                OC_value=0
            data=[[int(value),sp_unit_num],[OC_value,sp_out_unit_num]]
            # print(data)
            odd_ratio, p_value = stats.fisher_exact(data)
            # print(seq)
            seq_name=seq.split('_',2)[2]##suppose to be ful_intact_Cassandra_OS-LTR ==> Cassandra_OS-LTR
            status=seq.split('_',2)[1]##suppose to be ful_intact_Cassandra_OS-LTR ==> intact
            if p_value<0.05:
                dominate='dominate'
                line_dominate='\t'.join([sp,seq_name,status,str(value),str(sp_unit_num),str(OC_value),str(sp_out_unit_num),dominate])
            else:
                dominate='less'
            if dominate=='dominate':
                sp_group_unit_dom_out.write(line_dominate+'\n')
            line='\t'.join([sp,seq_name,status,str(value),str(sp_unit_num),str(OC_value),str(sp_out_unit_num),dominate])
            sp_group_unit_ful_out.write(line+'\n')
    sp_group_unit_dom_out.close()
    sp_group_unit_ful_out.close()
    print('sp_group done!')

    chr_group_unit_out=open(seq_stat_file.replace('01','02').split('.')[0]+'_dom_chr_long.txt','w')
    chr_group_unit_ful_out=open(seq_stat_file.replace('01','02').split('.')[0]+'_ful_chr_long.txt','w')
    chr_group_unit_out.write('Sp\tname\tstatue\tvalue\tunitnumber\tOCvalue\tOCunitvalue\tcondition'+'\n')
    chr_group_unit_ful_out.write('Sp\tname\tstatue\tvalue\tunitnumber\tOCvalue\tOCunitvalue\tcondition'+'\n')
    for index, row in chr_group_df.iterrows():
        sp=index
        sp_unit_num=int(chr_group_df.loc[sp,'unit_number'])
        sp_out_unit_num=int(OC_chr_group_df.loc[sp,'unit_number'])
        for column, value in row.items():
            if column=='unit_number':
                continue
            seq=column
            if seq in OC_chr_group_df:
                OC_value=int(OC_chr_group_df.loc[sp,seq])
            else:
                OC_value=0
            data=[[int(value),sp_unit_num],[OC_value,sp_out_unit_num]]
            # print(data)
            odd_ratio, p_value = stats.fisher_exact(data)
            # print(seq)
            seq_name=seq.split('_',2)[2]##suppose to be ful_intact_Cassandra_OS-LTR ==> Cassandra_OS-LTR
            status=seq.split('_',2)[1]##suppose to be ful_intact_Cassandra_OS-LTR ==> intact
            if p_value<0.05:
                dominate='dominate'
                line_dominate='\t'.join([sp,seq_name,status,str(value),str(sp_unit_num),str(OC_value),str(sp_out_unit_num),dominate])
            else:
                dominate='less'
            if dominate=='dominate':
                chr_group_unit_out.write(line_dominate+'\n')
            line='\t'.join([sp,seq_name,status,str(value),str(sp_unit_num),str(OC_value),str(sp_out_unit_num),dominate])
            chr_group_unit_ful_out.write(line+'\n')
    chr_group_unit_out.close()
    chr_group_unit_ful_out.close()
    print('chr_group done!')

if __name__ == "__main__":

    os.chdir(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult')

    # os.makedirs('Centro_region',exist_ok=True)
    # os.chdir(r'./Centro_region')
    seq_stat_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_centro_fullLTR_freq.csv"
    out_centro_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_OUTcentro_fullLTR_freq.csv"

    main(seq_stat_file,out_centro_file)

    seq_stat_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_centro_fullLTR_length.csv"
    out_centro_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_OUTcentro_fullLTR_length.csv"

    main(seq_stat_file,out_centro_file)