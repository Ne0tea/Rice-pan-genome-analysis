'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-25 20:09:19
LastEditors: Ne0tea
LastEditTime: 2024-02-03 14:44:48
'''
import os
import re

def split_clusters(clusters):
    new_clusters = [[0]]
    for i in range(1, len(clusters)):
        current_start, current_end = clusters[i][0],clusters[i][1]
        prev_start, prev_end = clusters[i-1][0],clusters[i-1][1]
        if current_start - prev_end > 500:
            new_clusters.append([i])
        else:
            new_clusters[-1].append(i)
    return new_clusters


def find_long_consecutive_duplicates_ranges_ignore_sd(lst, threshold=3):
    ranges = []
    start_index = None
    count = 1
    for i in range(len(lst) - 1):
        if lst[i] == lst[i + 1] :
            if count == 1 :
                start_index = i
            count += 1
        else:
            if count > threshold:
                ranges.append((start_index, i))
            count = 1
    # Check for the last element if it forms a consecutive range
    if count > threshold :
        ranges.append((start_index, len(lst) - 1))
    return ranges


def main(pan_gff_file_list,pan_centro_bed_file,Out_file):
    gff_list=[]
    Of=open(Out_file,'w')
    fulll=0
    ord_num=0
    with open(pan_gff_file_list,'r') as gf:
        for i in gf:
            line=i.strip()
            gff_list.append(line)
    pattern=re.compile(r'\"Motif:(.*)\"')
    cluster_number=0
    chr_centro_bed_dic={}
    with open(pan_centro_bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            if line[0] in chr_centro_bed_dic:
                chr_centro_bed_dic[line[0]][line[0].split('_')[1]]=[int(line[1]),int(line[2])]
            else:
                chr_centro_bed_dic[line[0]]={line[0].split('_')[1]:[int(line[1]),int(line[2])]}
    for file in gff_list:
        target_motif=[]
        target_motif_name_List=[]
        sp=os.path.basename(file).split('.')[0].replace('rmfup_replaced','')
        print(sp,'Start running')
        with open(file,'r') as gffile:
            count=0
            for i in gffile:
                if i.startswith('#'):
                    continue
                line=i.strip().split()
                '''
                length=int(line[4])-int(line[3])
                '''
                start=int(line[3])
                end=int(line[4])
                seq_name=pattern.findall(line[9])[0]
                chr_name=line[0]
                ord=line[6]
                if '(' in seq_name:
                    continue
                # if (chr_name=='Chr07' and sp=='CW06') or int(line[3]) > chr_centro_bed_dic[chr_name][sp][1] \
                #     or int(line[4]) < chr_centro_bed_dic[chr_name][sp][0] :
                #     continue
                target_motif.append({seq_name:[chr_name,start,end,seq_name,'0',ord]})
                target_motif_name_List.append(seq_name.split('_')[0])
                count+=1
        cluster=find_long_consecutive_duplicates_ranges_ignore_sd(target_motif_name_List)
        fulll+=len(cluster)
        for i in cluster:
            cur_cluster=[list(target_motif[index].keys())[0] for index in range(i[0],i[1]+1)]
            # if 'SZ-22_I' in cur_cluster:
            # print(cur_cluster.count('SZ-22_LTR'))
            cur_range_list=[list(target_motif[index].values())[0][1:3]for index in range(i[0],i[1]+1)]
            new_cluster=split_clusters(cur_range_list)
            # print(new_cluster)
            for unit_index in new_cluster:
                cluster_item=[list(target_motif[index+i[0]].keys())[0] for index in unit_index]
                if cluster_item.count('SZ-22_LTR') < 3:
                    continue
                cluster_number+=1
                # print(cluster_item)
                for inner_index in unit_index:
                    # print(inner_index)
                    ord_num+=1
                    cur_index=i[0]+inner_index
                    # print(list(target_motif[cur_index].values())[0])
                    line='\t'.join([str(unit) for unit in list(target_motif[cur_index].values())[0]])
                    # print(line)
                    unit_type=list(target_motif[cur_index].values())[0][3]
                    cluster_line=line+'\t'+'Cluster'+str(cluster_number)+'\t'+'Cluster'+str(cluster_number)+'_'+unit_type+'|'+str(ord_num)
                    # print(cluster_line)
                    Of.write(cluster_line+'\n')
        print(ord_num)
        print(sp,'done')
    # print(fulll)
    Of.close()


if __name__ == "__main__":
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff_rmfup.list"
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_periregion_70m.bed"
    Out_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22_cluster_distribution.bed"
    centro_state=0
    if int(centro_state) == 0:
        centro_stat=False
    else:
        centro_stat=True
    main(pan_gff_file_list,pan_centro_bed_file,Out_file)