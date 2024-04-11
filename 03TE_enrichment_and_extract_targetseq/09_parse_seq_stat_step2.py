'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-19 17:39:59
LastEditors: Ne0tea
LastEditTime: 2023-12-07 17:11:02
'''
import sys
import pandas as pd
import os
def seq_plus(x):
    final_dic={}
    for column_name, item in x.iteritems():
        if column_name == 'chr' or column_name=='sp':
            continue
        num_series=item.apply(lambda y:y.split('/')[0])
        num=num_series.astype(int).sum()
        length_series=item.apply(lambda y:y.split('/')[1])
        length=length_series.astype(int).sum()
        final_dic[column_name]=str(num)+'/'+str(length)
    final_series=pd.Series(final_dic).to_frame().T
    # final_series.index = pd.Series([i+'_'+sp])
    # print(final_series)
    return final_series

def main(seq_stat_file):
    seq_df=pd.read_csv(seq_stat_file)
    seq_df.columns = ['chr_sp'] + seq_df.columns[1:].tolist()

    seq_df[['chr','sp']]=seq_df['chr_sp'].str.split('_',expand=True)
    seq_df=seq_df.drop('chr_sp',axis=1)
    # sp_group_df=seq_df.drop('chr',axis=1)
    os.makedirs('parse_seq_step_tmp',exist_ok=True)
    sp_group_df=seq_df.groupby(['sp','chr']).apply(seq_plus)
    sp_group_df.to_csv(os.path.join('parse_seq_step_tmp','sp_chr_centro_group_ful.csv'))
    sp_group_df=pd.read_csv(os.path.join('parse_seq_step_tmp','sp_chr_centro_group_ful.csv'))
    sp_group_df=sp_group_df.loc[:,~sp_group_df.columns.str.contains("^Unnamed")]
    sp_group_df=sp_group_df.set_index('sp')
    os.makedirs('repeat_anno_group_by_chr',exist_ok=True)
    for i in sp_group_df['chr'].unique():
        c_df=sp_group_df[sp_group_df['chr']==i]
        # print(c_df)
        c_df=c_df.drop('chr',axis=1)
        freq_sp_df=pd.DataFrame()
        length_sp_df=pd.DataFrame()
        for column_name, item in c_df.iteritems():
            line=item.str.split('/',expand=True)[1].to_frame().astype('int')
            line.columns = [column_name]
            # print((line[column_name]==0).all())
            if (line[column_name]==0).all():
                continue
            else:
                length_sp_df=pd.concat([length_sp_df,line],axis=1)
            line=item.str.split('/',expand=True)[0].to_frame()
            line.columns = [column_name]
            freq_sp_df=pd.concat([freq_sp_df,line],axis=1)
        print(freq_sp_df)
        freq_sp_df.to_csv(os.path.join('repeat_anno_group_by_chr','freq_'+i+'_sp_df.csv'))
        length_sp_df.to_csv(os.path.join('repeat_anno_group_by_chr','length_'+i+'_sp_df.csv'))
    print()

if __name__ == "__main__":

    os.chdir(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult')
    
    os.makedirs('Centro_region',exist_ok=True)
    os.chdir(r'./Centro_region')
    seq_stat_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\Chr_seq_centro_region_intact.csv"
    
    # os.makedirs('Out_Centro_region',exist_ok=True)
    # os.chdir(r'./Out_Centro_region')
    # seq_stat_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\Chr_seq_out_centro_region_intact.csv"

    main(seq_stat_file)