'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-03-14 20:01:12
LastEditors: Ne0tea
LastEditTime: 2024-03-14 21:29:05
'''
import pandas as pd

def process_string(cell_value):
    # 使用正则表达式，去除第三个'_'之前的内容，并且去除字符串最后的'.'之后的内容
    # print(cell_value)
    processed_value = cell_value.split('_', 3)[-1].rsplit('.', 3)[0]
    return processed_value

bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_asm_70_CR_Msu7_clean.bed'
phy0=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_phylogenic.txt'

Sp_phy={}
Sp_index=[]
with open(phy0,'r') as phyf:
    for i in phyf:
        line=i.strip().split()
        Sp_index.append(line[0]+'_centro')
        Sp_phy[line[0]+'_centro']=line[1]

df=pd.read_table(bed_file,names=['Chr_sp','Start','End','ID','Strand'])

df['ID'] = df['ID'].apply(process_string)
df[['Chr', 'Sp']] = df['Chr_sp'].str.split('_', 1, expand=True)
df = df[df['Chr'] == 'Chr06']
print(df)

grouped_counts = df.groupby(['Sp', 'ID']).size().reset_index(name='Count')
pivot_table = grouped_counts.pivot_table(index='Sp', columns='ID', values='Count', fill_value=0)
pivot_table = pivot_table.reindex(Sp_index)
# print(pivot_table)
zero_counts = (pivot_table == 0).sum()
# 删除30个空值以上的列
columns_to_drop = zero_counts[zero_counts > 20].index
# print(columns_to_drop)
pivot_table.drop(columns=columns_to_drop, inplace=True)

pivot_table.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\CR_gene_anno.csv', index=True)

print(pivot_table)