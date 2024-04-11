'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-01-25 19:46:30
'''
import sys
import logging
import re

import numpy as np
from matplotlib import pyplot as plt

def find_intersecting_ranges(range_list, target_range):
    intersecting_ranges = []
    insert_loc=0
    for r in range_list:
        if (r[0] <= target_range[1] and r[1] >= target_range[0]) and \
            not (r[0] <= target_range[0] and r[1] >= target_range[1]):
            if r[0] <= target_range[0]:
                insert_loc=r[1]-target_range[0]
            else:
                insert_loc=r[0]-target_range[0]
            intersecting_ranges.append(r)
    return intersecting_ranges,insert_loc

def main(monomer_file,LTR_range_file):
    LTR_range_dic={}
    insert_loc_list=[]
    with open(LTR_range_file,'r') as lf:
        for i in lf:
            line=i.strip().split()
            if line[0] in LTR_range_dic:
                LTR_range_dic[line[0]].append((int(line[1]),int(line[2])))
            else:
                LTR_range_dic[line[0]]=[(int(line[1]),int(line[2]))]
    with open(monomer_file,'r') as mf:
        for i in mf:
            line=i.strip().split()
            if line[0] in LTR_range_dic:
                cur_range_list=LTR_range_dic[line[0]]
            else:
                continue
            cur_kmer_range=(int(line[1]),int(line[2]))
            intersect,insert_loc=find_intersecting_ranges(cur_range_list,cur_kmer_range)
            if intersect:
                # print(intersect,cur_kmer_range,insert_loc)
                insert_loc_list.append(insert_loc)
    print(len(insert_loc_list))

    plt.figure()
    plt.hist(insert_loc_list,edgecolor='black', alpha=0.7,bins=75,color='#bc8a57')

    # 添加标题和标签
    plt.title('Frequency Distribution Histogram')
    plt.xlim(xmin=0, xmax=160)  
    plt.xlabel('Values')
    plt.ylabel('Frequency')

    # 显示图形
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\LTR_insert_SAT_loc.pdf') 
if __name__ == "__main__":
    # unit_intact_anno_file=sys.argv[1]#SZ-22_unit_centro.bed 
    monomer_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\monomer_anno.txt'
    # unit_monomer_anno_file=sys.argv[2]#
    LTR_range_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_unit_centroV2.bed'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    concensus_length={'SZ-22_I':2912,
                      'SZ-22_LTR':792,
                      'CRM-LTR_OS':908,
                      'CRM-I_OS':5933,
                      'RIRE7_I':5899,
                      'RIRE7_LTR':858,
                      'GYPSY1-I_OS':7594,
                      'GYPSY1-LTR_OS':653,
                      'RETRO2A_I':10730,
                      'RETRO2A_LTR':1082,
                      'RETRO2_I':10690,
                      'RETRO2_LTR':1092,
                      'Gypsy-43_OS-I':1907,
                      'Gypsy-43_OS-LTR':1509,
                      'RIRE8A_I':6055,
                      'RIRE8A_LTR':2867}
    main(monomer_file,LTR_range_file)