'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-27 22:54:26
LastEditors: Ne0tea
LastEditTime: 2024-02-04 23:22:00
'''
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import seaborn as sns
import pandas as pd

def calculate_min_distance(data):
    only_cluster=0
    category_ranges = {}
    for entry in data:
        start, end, category = data[entry][0], data[entry][1], data[entry][2]
        if category not in category_ranges:
            category_ranges[category] = []
        category_ranges[category].append((start, end))
    min_distances = []
    for category, ranges in category_ranges.items():
        ranges.sort()  # 按起点排序
        min_distance = float('inf')
        if len(ranges) > 1:
            for i in range(1, len(ranges)):
                distance = ranges[i][0] - ranges[i-1][1]
                min_distance = min(min_distance, distance)
        else:
            only_cluster+=1
            continue
        min_distances.append(min_distance)
    # print(only_cluster)
    return min_distances

def calculate_pairwise_distance(data):
    only_cluster=0
    category_ranges = {}
    for entry in data:
        start, end, category = data[entry][0], data[entry][1], data[entry][2]
        if category not in category_ranges:
            category_ranges[category] = []
        category_ranges[category].append((start, end))
    # print(len(category_ranges['Chr01_13-65']))
    pairwise_distances = []
    for category, ranges in category_ranges.items():
        ranges.sort()  # 按起点排序
        if len(ranges) > 1:
            for x in range(0, len(ranges)):
                for y in range(x+1, len(ranges)):
                    distance = ranges[y][0] - ranges[x][0]
                    pairwise_distances.append(distance)
        else:
            only_cluster+=1
            continue
    # print(pairwise_distances)
    return pairwise_distances

def main(bed_file,bin_windows_file,peri_region_file,full_LTR_file,centro_region_file):
    Rention=[]
    cluster_centro_dic={}
    cluster_peri_dic={}
    cluster_dic={}
    cluster_LTR_number={}
    cluster_LTR_number_list=[]
    cluster_length={}
    cluster_length_list=[]
    unit_length_dic={'SZ-22_LTR':792,'SZ-22_I':2912}
    chr_cluster_LTR_number_dic={}
    chr_intact_number_dic={}
    chr_intact_number_dic_sp={}
    ord=0
    chr_bin_range={}
    chr_peri_range={}
    with open(peri_region_file,'r') as prf:
        for i in prf:
            line=i.strip().split()
            chr_peri_range[line[0]]=[int(line[1]),int(line[2])]
    chr_centro_range={}
    with open(centro_region_file,'r') as crf:
        for i in crf:
            line=i.strip().split()
            chr_centro_range[line[0]]=[int(line[1]),int(line[2])]
    with open(bin_windows_file,'r') as wf:
        for i in wf:
            line=i.strip().split()
            if line[0] in chr_bin_range:
                chr_bin_range[line[0]].append(int(line[1]))
            else:
                chr_bin_range[line[0]]=[int(line[1])]
        # print(len(chr_bin_range['Chr01_13-65']))
    chr_bin_range={key: sorted(value) for key, value in chr_bin_range.items()}
    intact_bin_item={}
    intact_s_e={}
    intact_bin_In,intact_bin_Out,cluster_In,cluster_Out=0,0,0,0
    intact={}
    cluster={}
    process_loc=[]
    with open(full_LTR_file,'r') as flf:
        for i in flf:
            line=i.strip().split()
            cur_s=int(line[1])
            cur_e=int(line[2])
            chr_sp_windows=chr_bin_range[line[0]]
            bin_site=sorted(chr_sp_windows.copy() + [cur_s]).index(cur_s)
            intact_name='_'.join(line[3].split('_')[:7])
            # print(intact_name)
            if intact_name not in process_loc:
                if 'Outperi_Cluster' in line[3]:
                    cluster_Out+=1
                    if 'Out' in cluster:
                        cluster['Out']+=1
                    else:
                        cluster['Out']=1
                elif 'Inperi_Cluster' in line[3]:
                    cluster_In+=1
                    if cur_s > chr_centro_range[line[0]][0] and cur_s < chr_centro_range[line[0]][1]:
                        if 'Centro' in cluster:
                            cluster['Centro']+=1
                        else:
                            cluster['Centro']=1
                    else:
                        if 'Pericentro' in cluster:
                            cluster['Pericentro']+=1
                        else:
                            cluster['Pericentro']=1
                elif 'Outperi_intact' in line[3]:
                    intact_bin_Out+=1
                    if 'Out' in intact:
                        intact['Out']+=1
                    else:
                        intact['Out']=1
                elif 'Inperi_intact' in line[3]:
                    intact_bin_In+=1
                    if cur_s > chr_centro_range[line[0]][0] and cur_s < chr_centro_range[line[0]][1]:
                        if 'Centro' in intact:
                            intact['Centro']+=1
                        else:
                            intact['Centro']=1
                    else:
                        if 'Pericentro' in intact:
                            intact['Pericentro']+=1
                        else:
                            intact['Pericentro']=1
            process_loc.append(intact_name)
            if line[3].split('_')[5].startswith('Cluster'):
                continue
            if str(bin_site) in intact_bin_item:
                intact_bin_item[str(bin_site)].append(intact_name)
            else:
                intact_bin_item[str(bin_site)]=[intact_name]
            if intact_name in intact_s_e:
                new_s=min(cur_s,intact_s_e[intact_name][0])
                new_e=max(cur_e,intact_s_e[intact_name][1])
                intact_s_e[intact_name]=[new_s,new_e]
            else:
                intact_s_e[intact_name]=[cur_s,cur_e]
    intact = dict(sorted(intact.items()))
    cluster = dict(sorted(cluster.items()))
    print(intact,cluster)
    # print(intact_bin_In,intact_bin_Out,cluster_In,cluster_Out)
    for i in chr_intact_number_dic_sp:
        nchr=i.split('_')[0]
        if nchr in chr_intact_number_dic:
            chr_intact_number_dic[nchr].append(chr_intact_number_dic_sp[i])
        else:
            chr_intact_number_dic[nchr]=[chr_intact_number_dic_sp[i]]

    cluster_bin_item={}
    cluster2chr={}

    with open(bed_file,'r') as bf:
        for i in bf:
            ord+=1
            line=i.strip().split()
            cluster2chr[line[6]]=line[0]
            cur_s=int(line[1])
            cur_e=int(line[2])
            Rention_ratio=(cur_e-cur_s)/unit_length_dic[line[3]] if (cur_e-cur_s)/unit_length_dic[line[3]] <=1 else 1
            Rention.append(Rention_ratio)
            chr=line[0].split('_')[0]
            chr_sp_windows=chr_bin_range[line[0]]
            bin_site=sorted(chr_sp_windows.copy() + [cur_s]).index(cur_s)
            if line[6] in cluster_dic:
                new_s=min(cur_s,cluster_dic[line[6]][0])
                new_e=max(cur_e,cluster_dic[line[6]][1])
                cluster_dic[line[6]]=[new_s,new_e,line[0]]
            else:
                cluster_dic[line[6]]=[cur_s,cur_e,line[0]]
            if str(bin_site) in cluster_bin_item:
                cluster_bin_item[str(bin_site)].append(line[6])
            else:
                cluster_bin_item[str(bin_site)]=[line[6]]
            if cur_s > chr_peri_range[line[0]][0] and cur_s < chr_peri_range[line[0]][1]:
                if line[6] in cluster_peri_dic:
                    new_s=min(cur_s,cluster_peri_dic[line[6]][0])
                    new_e=max(cur_e,cluster_peri_dic[line[6]][1])
                    cluster_peri_dic[line[6]]=[new_s,new_e,line[0]]
                else:
                    cluster_peri_dic[line[6]]=[cur_s,cur_e,line[0]]

            if cur_s > chr_centro_range[line[0]][0] and cur_s < chr_centro_range[line[0]][1]:
                if line[6] in cluster_centro_dic:
                    new_s=min(cur_s,cluster_centro_dic[line[6]][0])
                    new_e=max(cur_e,cluster_centro_dic[line[6]][1])
                    cluster_centro_dic[line[6]]=[new_s,new_e,line[0]]
                else:
                    cluster_centro_dic[line[6]]=[cur_s,cur_e,line[0]]
            if line[6] in cluster_LTR_number:
                if line[3] == 'SZ-22_LTR':
                    cluster_LTR_number[line[6]]+=1
            else:
                if line[3] == 'SZ-22_LTR':
                    cluster_LTR_number[line[6]]=1

    pairwise_dist=calculate_pairwise_distance(cluster_dic)
    pairwise_peri_dist=calculate_pairwise_distance(cluster_peri_dic)
    pairwise_centro_dist=calculate_pairwise_distance(cluster_centro_dic)

    cluster_pericentrolength_dic={}
    for i in cluster_dic:
        chr=cluster_dic[i][2].split('_')[0]
        c_unit_length=cluster_dic[i][1]-cluster_dic[i][0]
        cluster_length[i]=c_unit_length
        cluster_length_list.append(int(c_unit_length))
        if cluster_LTR_number[i] == 10:
            print(i)
        cluster_LTR_number_list.append(int(cluster_LTR_number[i]))
        if cluster_dic[i][2] in cluster_pericentrolength_dic:
            cluster_pericentrolength_dic[cluster_dic[i][2]].append(c_unit_length)
        else:
            cluster_pericentrolength_dic[cluster_dic[i][2]]=[c_unit_length]
        if chr in chr_cluster_LTR_number_dic:
            chr_cluster_LTR_number_dic[chr]+=int(cluster_LTR_number[i])
        else:
            chr_cluster_LTR_number_dic[chr]=int(cluster_LTR_number[i])
    cluster_bin_length={} #calculate median of cluster in bin
    cluster_bin_num={}
    for x in cluster_bin_item:
        bin_length_set=[]
        bin_num_set=[]
        # print(cluster_bin_item[x])
        for y in set(cluster_bin_item[x]):
            bin_length_set.append(cluster_length[y])
            bin_num_set.append(cluster_LTR_number[y])
        cluster_bin_length[x]=sum(bin_length_set)
        cluster_bin_num[x]=sum(bin_num_set)

    intact_bin_length={} #calculate median of cluster in bin
    intact_bin_num={}
    for x in intact_bin_item:
        bin_length_set=[]
        bin_num_set=0
        for y in set(intact_bin_item[x]):
            bin_length_set.append(intact_s_e[y][1]-intact_s_e[y][0])
            bin_num_set+=2
        intact_bin_length[x]=sum(bin_length_set)
        intact_bin_num[x]=bin_num_set

    for i in range(0,len(chr_bin_range['Chr01_13-65'])):
        if str(i) not in cluster_bin_length:
            cluster_bin_length[str(i)]=0
        if str(i) not in cluster_bin_num:
            cluster_bin_num[str(i)]=0
        if str(i) not in intact_bin_length:
            intact_bin_length[str(i)]=0
        if str(i) not in intact_bin_num:
            intact_bin_num[str(i)]=0
    cluster_length_dic={}
    cluster_LTR_number_list=np.array(cluster_LTR_number_list)
    for i in set(cluster_LTR_number_list):
        index=np.where(cluster_LTR_number_list==i)[0].tolist()
        # print(index)
        cluster_length_dic[i]=[cluster_length_list[x] for x in index]
    draw_data_list=[cluster_length_dic,cluster_length_dic,chr_cluster_LTR_number_dic,chr_intact_number_dic,
                    pairwise_peri_dist,pairwise_centro_dist,cluster_bin_num]
    title_list=['Cluster_LTR Vs Cluster_length','Cluster LTR porpotion','Chr Cluster LTR number','Intact and Cluster ratio',
                'Pairwise dist Peri','Pairwise dist Centro','Intact and Cluster number']
    fig = plt.figure(figsize=(12, 12))
    gs = GridSpec(3, 3)
    # draw total pattern plot
    for data, title in zip(draw_data_list, title_list):
        row_num=title_list.index(title)//3
        col_num=title_list.index(title)%3
        ax = plt.subplot(gs[row_num, col_num])
        if title=='Chr Cluster LTR number':
            x_values = list(data.keys())
            y_values = list(data.values())
            ax.scatter(x_values, y_values, label='散点图')
            ax.set_title('{} Distribution'.format(title))
        elif title=='Intact and Cluster ratio':
            data1=list(intact.values())
            data2=list(cluster.values())
            categories=['Centro','Out','Pericentro']
            ax.bar(categories, data1, label='Intact')
            ax.bar(categories, data2, bottom=data1, label='Cluster')
            ax.set_title('{}'.format(title))
            ax.legend()
        elif title=='Cluster LTR porpotion':
            datap=[len(data[x]) for x in data]
            # print(list(data.keys()))
            ax.pie(datap,labels=list(data.keys()))
            ax.legend()
            ax.axis('equal')
        elif title=='Intact and Cluster number':
            ax = plt.subplot(gs[2, :])
            sorted_data = dict(sorted(data.items(), key=lambda x: int(x[0])))
            categories = list(sorted_data.keys())
            values = list(sorted_data.values())
            p1=ax.scatter(categories, values,color="#bc6c25")
            p2=ax.plot(categories, values, linestyle='--', color='#bc6c25')

            ax2 = ax.twinx()
            sorted_data = dict(sorted(intact_bin_num.items(), key=lambda x: int(x[0])))
            categories = list(sorted_data.keys())
            values = list(sorted_data.values())
            p3=ax2.scatter(categories, values,color="#88ada6")
            p4=ax2.plot(categories, values, linestyle='--', color='#88ada6')
            p=ax.axvspan('80','160', facecolor='gray', alpha=0.3, label='Peri_Centromere')
            p_centro=ax.axvspan('105','120', facecolor='gray', alpha=0.5, label='Centromere')
            ax.set_title('{} Distribution'.format(title))
            average_interval = len(categories) // 10
            selected_labels = categories[::average_interval]
            ax.set_xticks(selected_labels)
            ax.legend([p1,p3,p,p_centro],['cluster LTR number in bin', 'intact LTR number in bin','Pericentromere region','Centromere region'])
        elif title=='Chr_Intact_number':
            df = pd.DataFrame([(k, v) for k, vs in data.items() for v in vs], columns=['Chr', 'Intact number'])
            sns.boxplot(x='Chr', y='Intact number', data=df, ax=ax, width=0.5,
                        order=['Chr01','Chr02','Chr03','Chr04','Chr05','Chr06',
                               'Chr07','Chr08','Chr09','Chr10','Chr11','Chr12'],showfliers=False)
            ax.set_title('{} Distribution'.format(title))
        elif title=='Cluster_LTR Vs Cluster_length':
            datap=list(cluster_length_dic.values())
            # print(datap)
            sns.violinplot(data=datap)
            fit = np.polyfit(cluster_LTR_number_list, cluster_length_list, 1)

            fit_x = np.linspace(min(cluster_LTR_number_list), max(cluster_LTR_number_list), 100)
            fit_y = np.polyval(fit, fit_x)
            # print(fit_x,fit_y)
            ax.plot(fit_x-3, fit_y, label='Trend line', color='#003049', linestyle='--')
            x_line = np.linspace(min(cluster_LTR_number_list), max(cluster_LTR_number_list), 100)
            y_line = 2248 * x_line
            ax.plot(x_line-3, y_line, color='#780000', linestyle='--', label='intact SZ-22 length')
            ax.legend()
            ax.set_title('{} Frequency Distribution'.format(title))
            ax.set_ylabel('Cluster length')
            ax.set_xticklabels(list(cluster_length_dic.keys()))
        else:
            ax.hist(data, bins=100, color='blue', alpha=0.7)
            ax.set_title('{} Frequency Distribution'.format(title))
            ax.set_ylabel('Frequency')
            x_ticks = np.linspace(min(data), max(data), 5)
            ax.set_xticks(x_ticks)

        if 'dist' in title and 'Pairwise_dist_I' != title:
            ax.set_xlabel('Distance')
            ax.set_xticklabels(['{:.4f}Mb'.format(i/1e6) for i in x_ticks])
        elif 'ratio' in title:
            ax.set_xlabel('Ratio')
        elif 'length' in title and 'Chr' not in title and 'loc' not in title and 'pericentro' not in title:
            ax.set_xlabel('Cluster LTR num')
        elif 'Pairwise_dist_I' == title:
            x_ticks = np.linspace(min(data), max(data), 5)
            ax.set_xlabel('Distance')
            # ax.set_xticklabels(['{:.4f}Mb'.format(i/1e6) for i in x_ticks])
        print(title,'plot finish')
        ax.grid(True)

    plt.tight_layout()  # 用于调整子图之间的间距，使得图形更美观
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\SZ-22_cluster_pattern.pdf') 


if __name__ == "__main__":
    cluster_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22_cluster_distribution.bed'
    bin_windows_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_sp.windows'
    centro_region_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    peri_region_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_periregion_70m.bed'
    full_LTR_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\merge_full_LTR.bed'
    main(cluster_anno_file,bin_windows_file,peri_region_file,full_LTR_file,centro_region_file)