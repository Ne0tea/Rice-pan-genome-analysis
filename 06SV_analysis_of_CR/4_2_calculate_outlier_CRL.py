'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-03-30 22:55:27
LastEditors: Ne0tea
LastEditTime: 2024-04-01 20:55:21
'''
import pandas as pd
from scipy import stats
import numpy as np

def outlier_evaluate(value,vlaue_list):
    mean_value = np.mean(vlaue_list)
    variance = np.var(vlaue_list)
    std_deviation = np.std(vlaue_list)
    value_to_measure = value
    deviation_from_mean = value_to_measure - mean_value
    z_score = deviation_from_mean / std_deviation
    return z_score

def identify_outliers(group):
    data=group.to_list()
    outliers=group.apply(outlier_evaluate, args=([data]))
    print(outliers)
    return outliers

def main(data,phy_file,outfile):
    wild_list=['obar','ogla']
    # geng_list=['or4','or3','trp','tmp']
    geng_list=['trp','tmp']
    aus_list=['aus']
    # xian_list=['XI2','XI3','XI1B','XI1A','or1','or2']
    xian_list=['XI2','XI3','XI1B','XI1A']
    phy_classfly={}
    with open(phy_file,'r') as phyf:
        for i in phyf:
            line=i.strip().split()
            if line[1] in wild_list:
                phy_classfly[line[0]]='Wild'
            elif line[1] in geng_list:
                phy_classfly[line[0]]='Geng'
            elif line[1] in aus_list:
                phy_classfly[line[0]]='Aus'
            elif line[1] in xian_list:
                phy_classfly[line[0]]='Xian'
    df = pd.read_table(data,names=['Chr_sp','S','E'])
    df['S'] = df['S'].astype(int)
    df['E'] = df['E'].astype(int)
    df['Length']=df['E']-df['S']
    df[['Chr', 'Sp']] = df['Chr_sp'].str.split('_', expand=True)
    df['Type'] = df['Sp'].map(phy_classfly)
    df=df[df['Type'].isin(['Geng','Xian'])]
    # 定义一个函数来标识离群值
    outliers = df.groupby(['Chr', 'Type'])['Length'].apply(identify_outliers)
    # print(outliers)
    df['isOut'] = outliers.astype(str)
    Out_df=df[df['isOut']!='False']
    Out_df=Out_df.reset_index(drop=True)

    Out_df.to_csv(outfile,index=False)
    print(Out_df)
    # print(Out_df[Out_df['isOut']=='Loss']['Sp'].value_counts())
    '''
    ooo=df[(df['Chr']=='Chr01') & (df['Type']=='Wild')]
    print(ooo)
    outliers = ooo.groupby(['Chr', 'Type'])['Length'].apply(identify_outliers)
    '''

if __name__ == "__main__":
    CR_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    phy_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_phylogenic.txt'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\Out_CRL.csv'
    main(CR_file,phy_file,outfile)