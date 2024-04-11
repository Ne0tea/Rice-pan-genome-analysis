'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-02-05 13:50:06
LastEditors: Ne0tea
LastEditTime: 2024-02-26 20:36:54
'''
import numpy as np
import re
import matplotlib.pyplot as plt
import pandas as pd
import PyComplexHeatmap as pch
from sklearn.preprocessing import MinMaxScaler

def parse_dist(file,order):
    line_number=0
    with open(file,'r') as df:
        ori_order=[]
        for i in df:
            line_number+=1
            if line_number <8:
                continue
            elif line_number==8:
                len_creat=len(set(i.strip().split()))
                continue
            line=i.strip().split()
            if 'data_matrix' in locals():
                # print(data_matrix.shape[1],len(order))
                cur_list=['nan']*(line_number-9)+list(line[:-2])+['nan']*(len_creat-len(line[:-2])-line_number+9)
                cur_arry=np.array([cur_list])
                # print(cur_arry)
                data_matrix=np.append(data_matrix,cur_arry,axis=0)
                ori_order.append(line[-2])
            else:
                # len_creat=len(order)
                data_matrix=np.array([['nan']*(line_number-9)+list(line[:-2])+['nan']*(len_creat-len(line[:-2])-line_number+9)])
                ori_order.append(line[-2])
    for i in order:
        if i not in ori_order:
            # print(i)
            ori_order.append(i)
            cur_list=['nan']*len_creat
            cur_arry=np.array([cur_list])
            data_matrix=np.append(data_matrix,cur_arry,axis=0)
            na_column = np.full((data_matrix.shape[0], 1), 'nan')
            data_matrix=np.append(data_matrix, na_column, axis=1)
    # np.savetxt(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\1111_LTR_matrix.txt', data_matrix, delimiter='\t',fmt = '%s')
    # print(data_matrix)

    # copy upper triangle
    upper_indices = np.triu_indices_from(data_matrix, k=1)
    data_matrix[upper_indices[1], upper_indices[0]] = data_matrix[upper_indices]
    # np.savetxt(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\2222_LTR_matrix.txt', data_matrix, delimiter='\t',fmt = '%s')

    # resort matrix according to give order
    # print(order)
    mapping = {value: index for index, value in enumerate(ori_order)}
    sorted_matrix = data_matrix[[mapping[i] for i in order]]
    sorted_matrix = sorted_matrix[:, [mapping[i] for i in order]]
    return sorted_matrix

def calculate_leaf_dist_array(x,genetree,gene_order):
    if x not in genetree:
        return [1]*len(gene_order)
    x_length=(genetree&x).dist
    cur_dist=[]
    # print(x)
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

def main(dist_dir,bed_file,outfile,centro_LTR_file=0):
    scaler = MinMaxScaler()
    data = pd.read_csv(bed_file, sep='\t', header=None)
    sorted_data = data.sort_values(by=[0, 1])
    fourth_column_order = sorted_data[3].tolist()
    centro_ltr_name=[]
    if centro_LTR_file:
        # start zoom in to centro region
        data = pd.read_csv(centro_LTR_file, sep='\t', header=None)
        sorted_data = data.sort_values(by=[0, 1])
        zoomIncentro = sorted_data[sorted_data[6]=='C'][3].tolist()
        # print(re.search(r'(.*)_', 'Chr09_NJ11_SZ-22_LTR_Inperi_intact_21_aw').group(1))
        fourth_column_order=[x for x in fourth_column_order if re.search(r'(.*)_', x).group(1)+'_1' in zoomIncentro]
        # print(fourth_column_order)
        outfile=outfile.split('.')[0]+'_centro.pdf'
    chr01_indices = [index for index, i in enumerate(fourth_column_order) if i.startswith('Chr01') ]
    chr01_txt=[i for i in fourth_column_order if i.startswith('Chr01') ]
    # print(chr01_txt)
    Chr_order = pd.DataFrame(fourth_column_order, columns=['Seq'])
    Chr_order['group'] = Chr_order['Seq'].str.split('_').str[0]
    Chr_order.index = fourth_column_order

    heatmapdata=parse_dist(dist_dir,fourth_column_order)
    heatmapdata = np.round(heatmapdata.astype(float),2)

    # max value of matrix to Na value
    Na_value=np.max(np.ma.masked_invalid(heatmapdata.astype(float)))
    heatmapdata[np.isnan(heatmapdata)] = Na_value
    heatmapdata[np.isinf(heatmapdata)] = Na_value

    # np.savetxt(outfile.split('.')[0]+'.matrix', heatmapdata, delimiter='\t', fmt='%.3f')
    chr01_sub_matrix = heatmapdata[np.ix_(chr01_indices, chr01_indices)]
    plt.figure(figsize=(10, 10))
    # np.savetxt(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_Chr01_fa_matrix.txt', chr01_sub_matrix, delimiter='\t', fmt='%.3f')
    np.savetxt(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_Chr01_LTR_matrix.txt', chr01_sub_matrix, delimiter='\t', fmt='%.3f')
    df_heatmap = pd.DataFrame(1-scaler.fit_transform(heatmapdata), columns=fourth_column_order)
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
                            # yticklabels_kws={'labelsize':7,'labelcolor':'blue'},
                            legend_kws={'vmax':1,'vmin':0.95},
                            legend_gap=5,cmap='coolwarm')
    plt.savefig(outfile, bbox_inches='tight')
    # plt.show()


if __name__ == '__main__':
    # tree_file=sys.argv[1]
    # centro_LTR_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\merge_full_LTR_LOCinfo.bed'
    
    numpy_dir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_fa500_pairwise_dist.npy'
    tree_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_SZ-22_LTR_extract_fa500.fasta.ali.trimal.tree"
    dist_dir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_SZ-22_LTR_extract_fa500.fasta.ali.trimal.distmat'
    bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_SZ-22_LTR_extract_fa500.bed'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\NJ11_SZ-22_LTR_extract_fa500.pdf'
    # numpy_dir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\SZ-22_pairwise_dist.npy'
    # dist_dir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\merge_NJ11_SZ-22_LTR.fasta.ali.trimal.distmat'
    # tree_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\merge_NJ11_SZ-22_LTR.fasta.ali.trimal.tree"
    # bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\merge_NJ11_SZ-22_LTR.bed'
    # outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_phylogeny\merge_NJ11_SZ-22_LTR.pdf'
    if 'centro_LTR_file' in locals():
        main(dist_dir,bed_file,outfile,centro_LTR_file)
    else:
        main(dist_dir,bed_file,outfile)

