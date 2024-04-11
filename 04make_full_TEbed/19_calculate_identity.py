'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-30 13:31:47
LastEditors: Ne0tea
LastEditTime: 2024-02-03 16:16:12
'''
from Bio.Align import PairwiseAligner
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.stats import mannwhitneyu
import seaborn as sns

def calculate_identity(seq1, seq2):
    aligner = PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    identity = best_alignment.score / max(len(seq1), len(seq2))
    return identity

def main(anno_file,fasta_file,bin_windows_file):
    cluster_set={}
    intact_set={}
    seq_dic={}
    LTR_list=[]
    cluster_LTR_number={}
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        seq_dic[seq.id]=seq.seq
    cluster_bin_item={}
    intact_bin_item={}
    chr_bin_range={}
    with open(bin_windows_file,'r') as wf:
        for i in wf:
            line=i.strip().split()
            if line[0] in chr_bin_range:
                chr_bin_range[line[0]].append(int(line[1]))
            else:
                chr_bin_range[line[0]]=[int(line[1])]
    chr_bin_range={key: sorted(value) for key, value in chr_bin_range.items()}
    with open(anno_file,'r') as af:
        for i in af:
            line=i.strip().split()
            chr_sp=line[0]
            start_site=int(line[1])
            chr_sp_windows=chr_bin_range[chr_sp]
            bin_site=sorted(chr_sp_windows.copy() + [start_site]).index(start_site)
            seq_name=line[3]
            loc_type=line[3].split('_')[5]
            if loc_type in cluster_LTR_number:
                if 'SZ-22_LTR' in line[3]:
                    cluster_LTR_number[loc_type]+=1
            else:
                if 'SZ-22_LTR' in line[3]:
                    cluster_LTR_number[loc_type]=1
            if loc_type.startswith('Cluster'):
                if loc_type in cluster_set:
                    cluster_set[loc_type].append(seq_name)
                else:
                    cluster_set[loc_type]=[seq_name]
                if str(bin_site) in cluster_bin_item:
                    cluster_bin_item[str(bin_site)].append(loc_type)
                else:
                    cluster_bin_item[str(bin_site)]=[loc_type]
            else:
                intact_num=line[0]+'_'+loc_type+str(line[3].split('_')[6])
                # if 'Outperi' in seq_name:
                #     continue
                if intact_num in intact_set:
                    intact_set[intact_num].append(seq_name)
                else:
                    intact_set[intact_num]=[seq_name]
                if str(bin_site) in intact_bin_item:
                    intact_bin_item[str(bin_site)].append(intact_num)
                else:
                    intact_bin_item[str(bin_site)]=[intact_num]
    # print(intact_bin_item)
    cluster_name_identity={}
    cluster_set_identity=[]
    cluster_set_identity1=[]
    cluster_set_identity2=[]
    cluster_set_identity3=[]
    for i in cluster_set:
        cur_cluster=cluster_set[i]
        cur_identity=[]
        for x in range(0,len(cur_cluster)):
            seq1=seq_dic[cur_cluster[x]]
            LTR_list.append(cur_cluster[x])
            for y in range(x,len(cur_cluster)):
                seq2=seq_dic[cur_cluster[y]]
                tmp_identity=calculate_identity(seq1,seq2)
                cur_identity.append(tmp_identity)
        if cluster_LTR_number[i] ==3:
            cluster_set_identity1.extend([np.median(cur_identity)]*len(cur_identity))
        elif cluster_LTR_number[i] == 4 :
            cluster_set_identity2.extend([np.median(cur_identity)]*len(cur_identity))
        elif cluster_LTR_number[i] >= 5 :
            cluster_set_identity3.extend([np.median(cur_identity)]*len(cur_identity))
        cluster_set_identity.extend([np.median(cur_identity)]*len(cur_identity))
        cluster_name_identity[i]=np.median(cur_identity)

    fail=0
    intact_name_identity={}
    intact_set_identity=[]
    for i in intact_set:
        if len(intact_set[i]) < 2:
            fail+=1
            continue
        intact_id1=intact_set[i][0]
        intact_id2=intact_set[i][1]
        LTR_list.append(intact_id1)
        LTR_list.append(intact_id2)
        intact_seq1=seq_dic[intact_id1]
        intact_seq2=seq_dic[intact_id2]
        cur_identity=calculate_identity(intact_seq1,intact_seq2)
        intact_set_identity.append(cur_identity)
        intact_name_identity[i]=cur_identity

    cluster_bin_identity={}
    for x in cluster_bin_item:
        bin_LTR_set=[]
        for y in cluster_bin_item[x]:
            bin_LTR_set.append(cluster_name_identity[y])
            # if x in cluster_bin_identity:
            #     cluster_bin_identity[x].append(cluster_name_identity[y])
            # else:
            #     cluster_bin_identity[x]=[cluster_name_identity[y]]
        cluster_bin_identity[x]=np.median(bin_LTR_set)

    intact_bin_identity={}
    for x in intact_bin_item:
        bin_LTR_set=[]
        for y in intact_bin_item[x]:
            if y not in intact_name_identity:
                continue
            bin_LTR_set.append(intact_name_identity[y])
            # if x in intact_bin_identity:
            #     intact_bin_identity[x].append(intact_name_identity[y])
            # else:
            #     intact_bin_identity[x]=[intact_name_identity[y]]
        intact_bin_identity[x]=np.mean(bin_LTR_set)
    heatmap_identity_LTR_list=[]
    for i in LTR_list:
        if 'Chr08_NJ11' in i:
            heatmap_identity_LTR_list.append(i)

    heatmap_identity_LTR_matrix = np.fromfunction(np.vectorize(lambda i, j: calculate_identity(seq_dic[heatmap_identity_LTR_list[int(i)]], seq_dic[heatmap_identity_LTR_list[int(j)]])), \
        (len(heatmap_identity_LTR_list), len(heatmap_identity_LTR_list)))
    row_clusters = sns.clustermap(heatmap_identity_LTR_matrix, method='average', cmap='viridis')
    row_order = row_clusters.dendrogram_row.reordered_ind
    col_order = row_clusters.dendrogram_col.reordered_ind
    # print(heatmap_identity_LTR_list)
    # print([heatmap_identity_LTR_list[i] for i in row_order])
    # 根据重新排序的索引对矩阵进行重新排列
    sorted_matrix = heatmap_identity_LTR_matrix[row_order][:, col_order]


    for i in range(0,len(chr_bin_range['Chr01_13-65'])):
        if str(i) not in cluster_bin_identity:
            cluster_bin_identity[str(i)]=0
        if str(i) not in intact_bin_identity:
            intact_bin_identity[str(i)]=0

    fig = plt.figure(figsize=(10, 8))
    gs = GridSpec(2, 2)
    # draw total pattern plot
    title_list=['Intact and Cluster LTR identity frequency','Intact and Cluster LTR identity','Bin Intact and Cluster LTR identity']
    for title in title_list:
        row_num=title_list.index(title)//4
        col_num=title_list.index(title)%2
        ax = plt.subplot(gs[row_num, col_num])
        if title=='Intact and Cluster LTR identity frequency':
            ax = plt.subplot(gs[0, 0])
            p1=ax.hist(intact_set_identity, bins=50, alpha=0.5, label='Intact', density=True)
            # ax.hist(cluster_set_identity, bins=50, alpha=0.5, label='Cluster', density=True)
            p2=ax.hist(cluster_set_identity1, bins=50, alpha=0.5, label='LTR 3', density=True)
            p3=ax.hist(cluster_set_identity2, bins=50, alpha=0.5, label='LTR 4', density=True)
            p4=ax.hist(cluster_set_identity3, bins=50, alpha=0.5, label='LTR o5', density=True)
            ax.set_title(title)
            ax.set_ylabel('Frequency')
            ax.legend()
            x_ticks = np.linspace(min(intact_set_identity), max(intact_set_identity), 5)
            ax.set_xticks(x_ticks)
            ax.set_xlim(0.4,1)
            ax.grid(True)
        # elif title=="Intact and Cluster LTR identity":
        #     ax = plt.subplot(gs[0, 1])

        #     ax.boxplot([intact_set_identity, cluster_set_identity], labels=['Intact', 'Cluster'])
        #     statistic, p_value = mannwhitneyu(intact_set_identity, cluster_set_identity, alternative='two-sided')
        #     if p_value < 0.05:
        #         ax.text(1.5, max(max(intact_set_identity), max(cluster_set_identity)), f'p = {p_value:.4f}\nSignificant', ha='center', va='center', color='red')
        #     else:
        #         ax.text(1.5, max(max(intact_set_identity), max(cluster_set_identity)), f'p = {p_value:.4f}\nNot Significant', ha='center', va='center', color='black')
        #     ax.set_ylabel('Identity')
        #     ax.set_title(title)
        #     ax.grid(True)

        elif title=='Bin Intact and Cluster LTR identity':
            ax = plt.subplot(gs[1, :])
            sorted_data = dict(sorted(cluster_bin_identity.items(), key=lambda x: int(x[0])))
            categories = list(sorted_data.keys())
            values = list(sorted_data.values())
            p1=ax.scatter(categories, values,color="#bc6c25")
            p2=ax.axvspan('80','160', facecolor='gray', alpha=0.3, label='Peri_Centromere')
            ax.plot(categories, values, linestyle='--', color='#bc6c25')

            sorted_data = dict(sorted(intact_bin_identity.items(), key=lambda x: int(x[0])))
            categories = list(sorted_data.keys())
            values = list(sorted_data.values())
            p3=ax.scatter(categories, values,color="#88ada6")
            ax.plot(categories, values, linestyle='--', color='#88ada6')
            average_interval = len(categories) // 10
            selected_labels = categories[::average_interval]
            ax.set_xticks(selected_labels)
            ax.set_title(title)
            ax.grid(False)
            ax.legend([p1,p3,p2],['Cluster LTR Identity', 'Intact LTR Identity','Pericentromere region'])
    plt.tight_layout()  # 用于调整子图之间的间距，使得图形更美观
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\SZ-22_cluster_Identity_pattern.pdf') 


    plt.figure(figsize=(12, 12))
    im=plt.imshow(sorted_matrix, cmap='RdBu_r', interpolation='nearest',vmin=0.2, vmax=1, aspect='auto')
    plt.xticks(np.arange(len(heatmap_identity_LTR_list)), [heatmap_identity_LTR_list[i] for i in row_order], rotation=90)
    plt.yticks(np.arange(len(heatmap_identity_LTR_list)), [heatmap_identity_LTR_list[i] for i in row_order])
    plt.colorbar(im, ax=ax, orientation='vertical')
    plt.title('Chr08_NJ11 LTR identity')
    plt.tight_layout()
    plt.savefig(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\SZ-22_cluster_Identity_pattern2.pdf') 
    
if __name__ == "__main__":
    full_LTR_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\merge_full_LTR.bed'
    seq_fasta_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\merge_full_LTR.fasta'
    bin_windows_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_sp.windows'
    main(full_LTR_file,seq_fasta_file,bin_windows_file)
