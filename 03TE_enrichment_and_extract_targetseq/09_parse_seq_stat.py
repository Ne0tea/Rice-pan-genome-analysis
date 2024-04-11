'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-18 18:06:31
LastEditors: Ne0tea
LastEditTime: 2023-12-01 21:19:08
'''
import sys
import re
import pandas as pd
import os

def seq_plus(x):
    final_dic={}
    # print(x)
    for column_name, item in x.iteritems():
        if column_name == 'chr' or column_name=='sp':
            continue
        num_series=item.apply(lambda y:y.split('/')[0])
        num=num_series.astype(int).sum()
        if column_name=='Cassandra_OS-LTR':
            print(num)
        length_series=item.apply(lambda y:y.split('/')[1])
        length=length_series.astype(int).sum()
        if column_name=='Cassandra_OS-LTR':
            print(str(num)+'/'+str(length))
        final_dic[column_name]=str(num)+'/'+str(length)
    final_series=pd.Series(final_dic).to_frame().T
    # final_series.index = pd.Series([i+'_'+sp])
    # print(final_series)
    return final_series

def main(seq_stat_file):
    os.makedirs('parse_seq_step_tmp',exist_ok=True)
    seq_df=pd.read_csv(seq_stat_file)
    # print(seq_df)
    seq_df.columns = ['chr_sp'] + seq_df.columns[1:].tolist()
    print(seq_df['chr_sp'])
    seq_df[['chr','sp']]=seq_df['chr_sp'].str.split('_',expand=True)
    seq_df=seq_df.drop('chr_sp',axis=1)

    chr_group_df=seq_df.drop('sp',axis=1)
    chr_group_df=chr_group_df.groupby('chr').apply(seq_plus)
    chr_group_df.to_csv(os.path.join('parse_seq_step_tmp','pan_chr_divide.csv'))

    # #test script
    # chr_order_chr_group_df=pd.DataFrame()
    # seq_order_chr_group_df=pd.DataFrame()
    # chr_group_df=pd.read_csv('E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_chr_divide.csv')
    # chr_group_df=chr_group_df.dropna()
    # chr_group_df=chr_group_df.loc[:,~chr_group_df.columns.str.contains("^Unnamed")]
    # chr_group_df.set_index('chr',inplace=True)
    # order_chr_group_df=chr_group_df

    # freq_chr_df=pd.DataFrame()
    # length_chr_df=pd.DataFrame()
    # for column_name, item in order_chr_group_df.iteritems():
    #     line=item.str.split('/',expand=True)[0].to_frame()
    #     line.columns = [column_name]
    #     freq_chr_df=pd.concat([freq_chr_df,line],axis=1)
    # for column_name, item in order_chr_group_df.iteritems():
    #     line=item.str.split('/',expand=True)[1].to_frame()
    #     line.columns = [column_name]
    #     length_chr_df=pd.concat([length_chr_df,line],axis=1)
    # colnames=sorted(freq_chr_df.columns.tolist(),reverse=True)
    # # print(olnames))
    # # print(freq_chr_df.shape[1])
    # freq_chr_df=freq_chr_df[colnames]
    # length_chr_df=length_chr_df[colnames]
    # freq_chr_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\freq_chr_df.csv')
    # length_chr_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\length_chr_df.csv')

    '''
    calucate rank of each seq
    for column_name, item in order_chr_group_df.iteritems():
        # print(item.str.split('/',expand=True)[0].map(int))
        line=item.str.split('/',expand=True)[0].map(int).rank(method='max',ascending=False).to_frame()
        line.columns = [column_name+'_rank']
        chr_order_chr_group_df=pd.concat([chr_order_chr_group_df,line],axis=1)
        # chr_order_chr_group_df[column_name+'_rank']=item.str.split('/',expand=True)[0].rank(ascending=False).map(int).to_frame().T
    chr_order_chr_group_df.loc['Col_sum'] = chr_order_chr_group_df.apply(lambda x: x.sum())
    chr_order_chr_group_df_T=chr_order_chr_group_df.T.sort_values('Col_sum')
    # print(chr_order_chr_group_df_T)
    for index, row in order_chr_group_df.iterrows():
        chr_num=index
        line=row.str.split('/',expand=True)[0].map(int).rank(method='max',ascending=False).to_frame().T
        line.index = [chr_num+'_rank']
        seq_order_chr_group_df=pd.concat([seq_order_chr_group_df,line])
        # seq_order_chr_group_df[chr_num+'_rank']=row.str.split('/',expand=True)[0].rank(ascending=False).map(int)
    # seq_order_chr_group_df['']
    seq_order_chr_group_df['Row_sum'] = seq_order_chr_group_df.apply(lambda x: x.sum(), axis=1)
    # print(seq_order_chr_group_df.T.sort_values('Col_sum'))
    # print(seq_order_chr_group_df)
    '''

    sp_group_df=seq_df.drop('chr',axis=1)
    sp_group_df=sp_group_df.groupby('sp').apply(seq_plus)
    sp_group_df.to_csv(os.path.join('parse_seq_step_tmp','pan_sp_divide.csv'))
    sp_group_df=pd.read_csv(os.path.join('parse_seq_step_tmp','pan_sp_divide.csv'))
    sp_group_df=sp_group_df.dropna()
    sp_group_df=sp_group_df.loc[:,~sp_group_df.columns.str.contains("^Unnamed")]
    # print(sp_group_df)
    sp_group_df.set_index('sp',inplace=True)
    order_sp_group_df=sp_group_df
    freq_sp_df=pd.DataFrame()
    length_sp_df=pd.DataFrame()
    for column_name, item in order_sp_group_df.iteritems():
        line=item.str.split('/',expand=True)[0].to_frame()
        line.columns = [column_name]
        freq_sp_df=pd.concat([freq_sp_df,line],axis=1)
    for column_name, item in order_sp_group_df.iteritems():
        line=item.str.split('/',expand=True)[1].to_frame()
        line.columns = [column_name]
        length_sp_df=pd.concat([length_sp_df,line],axis=1)
    colnames=sorted(freq_sp_df.columns.tolist(),reverse=True)
    freq_sp_df=freq_sp_df[colnames]
    length_sp_df=length_sp_df[colnames]
    freq_sp_df.to_csv(os.path.join('parse_seq_step_tmp','freq_sp_df.csv'))
    length_sp_df.to_csv(os.path.join('parse_seq_step_tmp','length_sp_df.csv'))
    # print(length_sp_df)


if __name__ == "__main__":
    #seq stat in every centromere generated by 08.py
    #example:chr_seq.stat
    os.chdir(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult')
    
    # os.makedirs('Centro_region',exist_ok=True)
    # os.chdir(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Centro_region')
    # seq_stat_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Chr_seq_centro_region_intact.csv"
    
    os.makedirs('Out_Centro_region',exist_ok=True)
    os.chdir(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Out_Centro_region')
    seq_stat_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Chr_seq_out_centro_region_intact.csv"

    # seq_stat_file=sys.argv[1]

    main(seq_stat_file)