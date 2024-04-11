'''
Descripttion: 
Out:获得specific unit之间的距离
Author: Ne0tea
version: 
Date: 2023-11-16 19:46:32
LastEditors: Ne0tea
LastEditTime: 2024-01-26 18:56:46
'''
import sys
import re
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import seaborn as sns

def distance_calulate(a,b):
    if a[0] > b[1]:
        return a[0]-b[1]
    elif a[1] < b[0]:
        return b[0]-a[1]
    else:
        print('There is the overlap ',a,'and',b)
        return 0
def find_nearest_distances(lst):
    nearest_distances = []

    for i in range(len(lst)):
        current_value = lst[i]
        remaining_values = lst[:i] + lst[i+1:]

        min_distance = min((distance_calulate(value,current_value)) for value in remaining_values)

        nearest_distances.append(min_distance)

    return nearest_distances

def main(unit_loc_file,spec_seq):
    loc_dic={}

    with open(unit_loc_file,'r') as loc_file:
        for i in loc_file:
            line=i.strip().split()
            if spec_seq not in line[3]:
                continue
            sp=line[0].split('_')[1]
            chr=line[0].split('_')[0]
            start=int(line[1])
            end=int(line[2])
            # location=line[5]
            if sp not in loc_dic:
                loc_dic[sp]={chr:[[start,end]]}
            else:
                if chr not in loc_dic[sp]:
                    loc_dic[sp][chr]=[[start,end]]
                else:
                    loc_dic[sp][chr].append([start,end])
    ful_seq_list=[]
    sns.set(style="whitegrid")  # 设置风格，这里选择了白色网格背景
    plt.figure(figsize=(60, 56))
    sp_number=1
    # print(loc_dic)
    for sp in loc_dic:
        sp_seq_list=[]
        plt.figure(figsize=(18, 10))
        chr_number=1
        for chr in loc_dic[sp]:
            cur_loc=loc_dic[sp][chr]
            if len(cur_loc)>1:
                pw_dl=find_nearest_distances(cur_loc)
            else:
                continue
            ful_seq_list.extend(pw_dl)
            sp_seq_list.extend(pw_dl)
            # # 绘制密度图
            # print(chr_number)
            plt.subplot(3, 4, chr_number)
            plt.hist(pw_dl, fill=True)
            plt.title('Density Plot')  # 设置标题
            plt.xlabel('Values')  # 设置 x 轴标签
            plt.ylabel('Density')  # 设置 y 轴标签
            data_count = len(pw_dl)
            text_to_display = f'Data Count: {data_count}'
            x_max, y_max = plt.gca().get_xlim()[1], plt.gca().get_ylim()[1]
            # text_x = x_max  # 在 x 轴上设置一点偏移量
            # x_max=2.1e6
            text_y = y_max * 0.9   # 在 y 轴上设置一点偏移量
            print(text_to_display)
            plt.text(x_max, text_y, text_to_display, fontsize=12, color='black', ha='right')
            # plt.xlim(0, 2.1e6)
            chr_number+=1
            # plt.show()
        plt.savefig(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\Ful_'+sp+r'_plots.pdf', format='pdf')
        plt.close()
        plt.subplot(8, 9, sp_number)
        plt.hist(sp_seq_list, fill=True)
        plt.title('Density Plot')  # 设置标题
        plt.xlabel('Values')  # 设置 x 轴标签
        plt.ylabel('Density')  # 设置 y 轴标签
        data_count = len(sp_seq_list)
        text_to_display = f'Data Count: {data_count}'
        x_max, y_max = plt.gca().get_xlim()[1], plt.gca().get_ylim()[1]
        # text_x = x_max  # 在 x 轴上设置一点偏移量
        text_y = y_max * 0.95   # 在 y 轴上设置一点偏移量
        
        # print(x_max, text_y)
        print(sp,text_to_display)
        # x_max=2.1e6
        plt.text(x_max, text_y*0.95, text_to_display, fontsize=9, color='black', ha='right')
        plt.text(x_max, text_y, sp, fontsize=9, color='black', ha='right')
        # plt.xlim(0, 2.1e6)
        sp_number+=1
    plt.tight_layout()
    plt.savefig('Ful_sp_plots.pdf', format='pdf')
    # plt.show()
    
    sns.set(style="whitegrid")  # 设置风格，这里选择了白色网格背景
    plt.figure(figsize=(8, 6))  # 设置图形大小

    # 绘制密度图
    plt.hist(ful_seq_list, fill=True)
    # plt.hist(ful_seq_list, bins=30, alpha=0.7, color='blue')
    plt.title('Density Plot')  # 设置标题
    plt.xlabel('Values')  # 设置 x 轴标签
    plt.ylabel('Density')  # 设置 y 轴标签
    data_count = len(ful_seq_list)
    text_to_display = f'Data Count: {data_count}'
    x_max, y_max = plt.gca().get_xlim()[1], plt.gca().get_ylim()[1]
    # text_x = x_max  # 在 x 轴上设置一点偏移量
    text_y = y_max * 0.95   # 在 y 轴上设置一点偏移量

    print(text_to_display)
    plt.text(x_max, text_y, text_to_display, fontsize=12, color='black', ha='right')
    plt.savefig('Ful_density.pdf', format='pdf')

if __name__ == "__main__":
    # unit_loc_file=sys.argv[1]
    # sp,chr,ltr_name,str(unit_start),str(unit_end),'centro'
    unit_loc_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_unit_centroV2.bed'
    # spec_seq=sys.argv[2]
    spec_seq='SZ-22'

    main(unit_loc_file,spec_seq)