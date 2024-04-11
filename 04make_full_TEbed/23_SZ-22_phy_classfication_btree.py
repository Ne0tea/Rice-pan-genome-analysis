'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-02-05 13:50:06
LastEditors: Ne0tea
LastEditTime: 2024-02-18 15:30:36
'''
import numpy as np
from ete3 import Tree
import os
import matplotlib.pyplot as plt
import pandas as pd
import PyComplexHeatmap as pch
from sklearn.preprocessing import MinMaxScaler
def calculate_leaf_dist_array(x,genetree,gene_order):
    if x not in genetree:
        return [1]*len(gene_order)
    x_length=(genetree&x).dist
    cur_dist=[]
    print(x)
    for y in gene_order:
        if y in genetree:
            y_length=(genetree&y).dist
            # print(y)
            if x==y:
                cur_dist.append(0)
            else:
                cur_dist.append((genetree&x).get_distance(y)-x_length-y_length)
        else:
            cur_dist.append(1)
    return cur_dist


if __name__ == '__main__':
    # tree_file=sys.argv[1]
    # numpy_dir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_fa500_pairwise_dist.npy'
    # tree_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_SZ-22_LTR_extract_fa500.fasta.ali.trimal.tree"
    # bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_SZ-22_LTR_extract_fa500.bed'
    numpy_dir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\SZ-22_pairwise_dist.npy'
    tree_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\merge_NJ11_SZ-22_LTR.fasta.ali.trimal.tree"
    bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\merge_NJ11_SZ-22_LTR.bed'
    if not os.path.exists(numpy_dir):
        gene_tree=Tree(tree_file)
        gene_name=gene_tree.get_leaf_names()
        # 读取txt文件，假设文件以空格分隔，没有列名
        data = pd.read_csv(bed_file, sep='\t', header=None)

        # 按照第一列和第二列的值排序
        sorted_data = data.sort_values(by=[0, 1])

        # 返回第四列的顺序
        fourth_column_order = sorted_data[3].tolist()
        data_np=np.empty((0,len(fourth_column_order)),float)
        for i in fourth_column_order:
            # print(i)
            cur_dis_array=calculate_leaf_dist_array(i,gene_tree,fourth_column_order)
            cur_arry=np.array([cur_dis_array])
            data_np=np.append(data_np,cur_arry,axis=0)
        np.save('SZ-22_pairwise_dist.npy',data_np)
        print("Array was created.")
    else:
        scaler = MinMaxScaler()
        ld=np.load(numpy_dir)
        # print(ld)
        ld[np.tril_indices_from(ld)] = 0
        data = pd.read_csv(bed_file, sep='\t', header=None)
        sorted_data = data.sort_values(by=[0, 1])
        fourth_column_order = sorted_data[3].tolist()
        chr01_indices = [i for i, index in enumerate(fourth_column_order)  ]
        chr01_sub_matrix = ld[np.ix_(chr01_indices, chr01_indices)]
        Chr_order = pd.DataFrame(fourth_column_order, columns=['Seq'])
        Chr_order['group'] = Chr_order['Seq'].str.split('_').str[0]
        Chr_order.index = fourth_column_order
        # np.savetxt(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_fa_matrix.txt', chr01_sub_matrix, delimiter='\t')
        # np.savetxt(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_LTR_matrix.txt', chr01_sub_matrix, delimiter='\t')
        # 绘制热图
        plt.figure(figsize=(12, 12))
        df_heatmap = pd.DataFrame(scaler.fit_transform(ld), columns=fourth_column_order,dtype=float)
        df_heatmap.index = fourth_column_order
        cm = pch.ClusterMapPlotter(data=df_heatmap,
                                col_cluster=False,row_cluster=False,
                                label='genetic distance',
                                show_rownames=True,show_colnames=True,row_names_side='left',
                                col_split_gap=1,row_split_gap=1,
                                row_split_order=sorted(set(Chr_order.group.values)),
                                row_split=Chr_order.group,
                                col_split=Chr_order.group,
                                col_split_order=sorted(set(Chr_order.group.values)),
                                row_dendrogram=False,col_dendrogram=False,
                                # col_names_side= "top",
                                legend_gap=5,cmap='turbo',
                                yticklabels_kws={'labelsize':7,'labelcolor':'blue'})
        # plt.savefig("example0.pdf", bbox_inches='tight')
        plt.show()
