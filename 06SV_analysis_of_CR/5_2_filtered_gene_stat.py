'''
Descripttion: gene_file:gene注释的bed文件
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-04-03 10:32:50
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
pd.set_option('display.max_rows', None)
def draw_plot(plot_data):
    p1_data=plot_data[0]
    p2_data=plot_data[1]
    p3_data=plot_data[2]
    p4_data=plot_data[3]
    p5_data=plot_data[4]
    p6_data=plot_data[5]
    
    fig, axs = plt.subplots(2, 3, figsize=(10, 6))
    # 绘制直方图
    axs[0, 0].hist(p1_data, bins=50, density=True, color='skyblue')
    axs[0, 0].text(0.95, 0.95, f'n={len(p1_data)}', transform=axs[0, 0].transAxes, ha='right', va='top')
    axs[0, 0].set_title('Gene Number')
    # 绘制直方图
    axs[1, 0].hist(p2_data, bins=50, density=True, color='skyblue')
    axs[1, 0].set_title('Gene Length')

    # 绘制箱型图
    axs[0, 1].bar(p3_data[0],p3_data[1])
    axs[0, 1].set_xticklabels(p3_data[0],rotation=45)
    axs[0, 1].set_title('Gene number / Chr')

    # 绘制第二个箱型图
    axs[0, 2].boxplot(p4_data[1],labels=p4_data[0])
    axs[0, 2].set_xticklabels(p4_data[0],rotation=45)
    axs[0, 2].set_title('Gene length / Chr')

    # 绘制箱型图
    axs[1, 1].bar(p5_data[0],p5_data[1])
    axs[1, 1].set_title('Gene number / Phy')

    # 绘制第二个箱型图
    axs[1, 2].boxplot(p6_data[1],labels=p6_data[0])
    # axs[1, 1].xticks(range(len(p5_data[1])), p5_data[0])
    axs[1, 2].set_title('Gene length / Phy')
    # 调整布局
    plt.tight_layout()

    # 显示图形
    plt.show()

def get_len_list(group):
    return group['len'].tolist()
def read_file_info(gff3_file,phy_ord_file):
    phy_dic={}
    with open(phy_ord_file,'r') as gfile:
        for i in gfile:
            line=i.strip().split()
            phy_dic[line[0]]=line[2]

    draw_data=[]
    df=pd.read_table(gff3_file,names=['Chr_sp','source','type','start','end','null1','strand','null2','attributes'])
    df=df[df['type']=='gene']
    df['phy']=df['Chr_sp'].map(lambda x:phy_dic[x.split('_')[1]])
    df['len']=df['end']-df['start']
    df['Chr']=df['Chr_sp'].map(lambda x:x.split('_')[0])
    df['Sp']=df['Chr_sp'].map(lambda x:x.split('_')[1])
    gene_length = df['len'].to_list()
    gene_number = df.groupby('Sp').size()
    gene_number_value=gene_number.to_list()

    max_gene_SpChr,max_gene_num=df.groupby(['Sp','Chr']).size().idxmax(),df.groupby(['Sp','Chr']).size().max()
    print('max gene number insert in CR',max_gene_SpChr,max_gene_num)
    max_gene_phySp,max_gene_num=df.groupby(['phy','Sp']).size().idxmax(),df.groupby(['phy','Sp']).size().max()
    print('max gene number insert in Phy',max_gene_phySp,max_gene_num)
    phy_sp_maxg=df.groupby(['phy','Sp']).size().reset_index().groupby(['phy']).max().reset_index()
    phy_sp_maxg=phy_sp_maxg.set_index(phy_sp_maxg['phy'])
    grouped = df.groupby(['phy', 'Sp']).size()
    max_sp_index = grouped.groupby(level='phy').idxmax()
    phy_sp_maxg['value']=max_sp_index
    phy_sp_maxg = phy_sp_maxg.drop(columns=['phy','Sp'])
    phy_sp_maxg=phy_sp_maxg.rename(columns={ 0: 'MaxG', 'value': 'MaxSp'})
    print(phy_sp_maxg.reset_index())


    phy_gnum = df.groupby('phy').size()
    phy_sp_counts = df.groupby('phy')['Sp'].nunique()
    phy_gnum=phy_gnum/phy_sp_counts

    phy_glen_index = df.groupby('phy').apply(get_len_list).index.to_list()
    phy_glen = df.groupby('phy').apply(get_len_list)
    Chr_gnum = df.groupby('Chr').size()
    Chr_glen_index = df.groupby('Chr').apply(get_len_list).index.to_list()
    Chr_glen = df.groupby('Chr').apply(get_len_list)

    Chr_gnum_index = Chr_gnum.index.to_list()
    Chr_gnum_value = Chr_gnum.to_list()
    Chr_glen_matrix = np.array(Chr_glen)

    phy_gnum_index=phy_gnum.index.to_list()
    phy_gnum_value=phy_gnum.to_list()
    phy_glen_matrix=np.array(phy_glen)

    draw_data.append(gene_number_value)
    draw_data.append(gene_length)
    draw_data.append([Chr_gnum_index,Chr_gnum_value])
    draw_data.append([Chr_glen_index,Chr_glen_matrix])
    draw_data.append([phy_gnum_index,phy_gnum_value])
    draw_data.append([phy_glen_index,phy_glen_matrix])

    return draw_data

def main(gff3_file,phy_ord_file):
    draw_data=read_file_info(gff3_file,phy_ord_file)
    # draw_plot(draw_data)
    # print(draw_data)

if __name__ == "__main__":
    gff3_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_asm_70_CR_Msu7_clean_filtered2.gff3'
    phy_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_phylogenic.txt'
    main(gff3_file,phy_file)
