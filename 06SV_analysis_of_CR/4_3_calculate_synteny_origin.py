'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-03-31 15:11:14
LastEditors: Ne0tea
LastEditTime: 2024-04-06 19:52:11
'''
import pandas as pd
import numpy as np
from scipy.spatial import distance
from collections import Counter
import os
pd.set_option('display.max_rows', None)

def calculate_ratio(row):
    total = row.sum()
    return row / total

def scale_matrix_to_range(matrix, new_min, new_max):
    min_value = np.min(matrix)
    max_value = np.max(matrix)
    scaled_matrix = ((matrix - min_value) / (max_value - min_value)) * (new_max - new_min) + new_min
    return scaled_matrix

def closest_row_value(lst,df,phy):
    # 计算欧氏距离 np.sqrt计算的结果和distance结果一致
    distances = np.sqrt(((df.values - lst) ** 2).sum(axis=1))
    # distances2 = np.array([distance.euclidean(row, lst) for row in df.values])
    closest_index = df.index.to_list().index(phy)
    return distances[closest_index]
    # return df.index[closest_row_index]
    
def closest_row_index(lst,df,cur_phy,phy_dic):
    '''
    Geng and Xian bias according to sample bias in wild
    '''
    distances = np.array([distance.euclidean(row, lst) for row in df.values])
    # # 找到距离最小的行的索引
    row_index = df.index[distances.argsort()[:5]]
    closest_phy = [ phy_dic[x] for x in row_index]
    # print(phy_dic[cur_phy],closest_phy)
    if closest_phy.count(phy_dic[cur_phy]) > 2:
        return 'Same_phy'
    else:
        return 'Vari_phy'

def main(data,phy_file,outdir):
    total_synteny_ori_df=pd.DataFrame()
    phydf=pd.read_table(phy_file,names=['Sp','ori_Phy'],index_col=0)

    df = pd.read_table(data,names=['Chr_sp','S','E','Strand','ID','family','Sp','Phy'])
    df[['Chr', 'Sp']] = df['Chr_sp'].str.split('_', expand=True)
    phy_dic={}
    for index,Sp in df[['Phy', 'Sp']].iterrows():
        if Sp.iloc[1] in phy_dic:
            continue
        phy_dic[Sp.iloc[1]]=Sp.iloc[0]
    tmp_family_dataset = df.groupby(['Chr', 'Phy','Sp'])['family'].value_counts().unstack(fill_value=0)
    family_monomer_matrix = tmp_family_dataset.groupby(level=['Chr', 'Phy']).median()
    family_monomer_matrix=family_monomer_matrix.reset_index()
    # print(family_monomer_matrix)
    ''' calculate average value of every phy
    family_count=df.groupby(['Chr', 'Phy'])['family'].value_counts()
    family_Sp_count = df.groupby(['Chr', 'Phy'])['Sp'].nunique()
    family_monomer_matrix=family_count/family_Sp_count
    '''
    for chr in family_monomer_matrix['Chr'].unique():
        cur_family_dataset_file_dir=os.path.join(outdir,chr+'_family_dataset.csv')
        cur_Sp_dataset_file_dir=os.path.join(outdir,chr+'_Sp_dataset.csv')
        cur_family_monomer=family_monomer_matrix[family_monomer_matrix['Chr']==chr]

        cur_family_monomer.set_index('Phy',inplace=True)
        cur_family_monomer.drop(columns=['Chr'], inplace=True)
        cur_family_monomer = cur_family_monomer.loc[:, (cur_family_monomer > 50).any()]
        cur_target_family=cur_family_monomer.columns.to_list()
        cur_family_monomer_dataset=cur_family_monomer

        '''convert to ratio'''
        cur_family_monomer_dataset=cur_family_monomer_dataset.apply(calculate_ratio, axis=1)
        print(cur_family_monomer_dataset)

        # cur_family_monomer_dataset.to_csv(cur_family_dataset_file_dir)
        cur_Sp_df=df[df['Chr']==chr]
        cur_Sp_monomer=cur_Sp_df.groupby(['Chr', 'Sp'])['family'].value_counts()
        cur_Sp_monomer.name='value'
        # print(cur_Sp_monomer)
        cur_Sp_monomer_df=pd.DataFrame(cur_Sp_monomer).reset_index()
        cur_Sp_monomer_df=cur_Sp_monomer_df[cur_Sp_monomer_df['family'].isin(cur_target_family)]
        # print(cur_Sp_monomer_df)
        cur_Sp_dataset=cur_Sp_monomer_df.pivot_table(index='Sp', columns='family', values='value', aggfunc='first',fill_value=0)

        '''convert to ratio'''
        cur_Sp_dataset=cur_Sp_dataset.apply(calculate_ratio, axis=1)
        # print(cur_Sp_dataset)

        cur_Sp_phy_monomer_vari=cur_Sp_dataset.apply(lambda row: closest_row_value(row.tolist(), cur_family_monomer_dataset,phy_dic[row.name]), axis=1)
        cur_Sp_closest_phy=cur_Sp_dataset.apply(lambda row: closest_row_index(row.tolist(), cur_Sp_dataset,row.name,phy_dic), axis=1)
        cur_Sp_dataset['monomer_Vari']=cur_Sp_phy_monomer_vari
        cur_Sp_dataset['monomer_Status']=cur_Sp_closest_phy
        cur_Sp_dataset['ori_Phy']=phydf['ori_Phy']
        cur_Sp_dataset['Chr_Sp']=chr+'_'+cur_Sp_dataset.index
        cur_Sp_dataset.set_index(cur_Sp_dataset['Chr_Sp'], inplace=True)
        # cur_Sp_dataset.to_csv(cur_Sp_dataset_file_dir)
        # total_synteny_ori_df=pd.concat([total_synteny_ori_df,cur_Sp_dataset[['monomer_Vari','monomer_Status','ori_Phy']]])
        print(cur_Sp_dataset)
    # total_synteny_ori_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\70_CR_synteny_ori.csv')
if __name__ == "__main__":
    CR_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\monomer_test.txt'
    phy_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_phylogenic.txt'
    outdir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR'
    main(CR_file,phy_file,outdir)