'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-23 14:55:08
LastEditors: Ne0tea
LastEditTime: 2024-01-24 22:12:22
'''
import sys
import logging
import re

import numpy as np
from matplotlib import pyplot as plt
 

def cal_distance(site,range1):
    ranges=min(range1)
    rangee=max(range1)
    if site > ranges and site < rangee:
        return 0
    elif site < ranges :
        return ranges-site
    elif site > rangee :
        return site-rangee

def find_in_LTR(range1, dic1):
    dic2={x:cal_distance(dic1[x],range1) for x in dic1}
    min_keys = sorted(dic2, key=dic2.get)[:2]
    min_type=[x.split('|')[0] for x in min_keys]
    min_statue=[x.split('|')[2] for x in min_keys]
    if min_type[0]==min_type[1]:
        loc='LTR' if 'LTR' in min_type[0] else 'Inner'
        if 'fulllength' in min_statue:
            # print('artigue')
            return loc+' | In LTR ful length'
        else:
            return loc+' | In LTR with loss'
    else:
        if 'loss' in min_statue:
            return 'Between LTR with loss'
        else:
            return 'Between LTR ful length'

def main(intact_file,monomer_file):
    global concensus_length
    pattern=re.compile(r'\"Motif:(.*)\"')
    LTR_unit_dic={}
    with open(intact_file,'r') as iff:
        count=0
        for i in iff:
            count+=1
            if i.startswith('#'):
                continue
            if 'RCH2' in i or 'SAT' in i or '(T)n' in i or 'LTR-25_OS' in i \
                or 'Cassandra_OS' in i or '(AT)n' in i:
                continue
            line=i.strip().split()
            LTR_name=line[3]
            start=int(line[1])
            end=int(line[2])
            if pattern.findall(line[11])[0] in concensus_length:
                cl=concensus_length[pattern.findall(line[11])[0]]
            else:
                cl=10
            LTR_statue='fulllength' if (end-start)/ cl > 0.8 else 'loss'
            Ltr_type=pattern.findall(line[11])[0]+'|'+str(count)+'|'+LTR_statue+'|'+str(start)
            if LTR_name in LTR_unit_dic:
                LTR_unit_dic[LTR_name][Ltr_type]=[start,end]
            else:
                LTR_unit_dic[LTR_name]={Ltr_type:[start,end]}
    processed_unit=[]
    tmp=[]
    with open(monomer_file,'r') as mff:
        inter_count=0
        LTR_count=0
        Inner_count=0
        between_count=0
        for i in mff:
            if i.startswith('#'):
                continue
            line=i.strip().split()
            LTR_name=line[3]
            start=int(line[8])
            end=int(line[9])
            if LTR_name in processed_unit:
                continue
            else:
                processed_unit.append(LTR_name)
            intersect_statue=False
            cur_unit_dic=LTR_unit_dic[LTR_name]
            cur_unit_start={key:cur_unit_dic[key][0] for key in cur_unit_dic}
            insert_statue=find_in_LTR([start,end],cur_unit_start)
            print(insert_statue)
            if 'In LTR' in insert_statue:
                inter_count+=1
                SAT_insert_site=start-int(line[1])
                tmp.append(SAT_insert_site)
                if insert_statue.startswith('LTR'):
                    LTR_count+=1
                elif insert_statue.startswith('Inner'):
                    Inner_count+=1
            else:
                between_count+=1
                SAT_insert_site=start-int(line[1])
                tmp.append(SAT_insert_site)
        
        plt.figure()
        plt.pie([LTR_count,Inner_count,between_count],labels=["LTR","Inner","Between"],colors=['#bc8a57','#b7c7a0','#88ada6'])
        plt.show()
        
        plt.figure()
        plt.hist(tmp,edgecolor='black', alpha=0.7,bins=100,color='#505050')

        # 添加标题和标签
        plt.title('Frequency Distribution Histogram')
        plt.xlim(xmin=0, xmax=4500)  
        plt.xlabel('Values')
        plt.ylabel('Frequency')

        # 显示图形
        plt.show()
        logging.info(f"SAT insert in LTR:{inter_count}")
        logging.info(f"SAT insert between LTR:{between_count}")
if __name__ == "__main__":
    # unit_intact_anno_file=sys.argv[1]#SZ-22_unit_centro.bed 
    unit_intact_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_unit_centroV2.detailed.bed'
    # unit_monomer_anno_file=sys.argv[2]#
    unit_monomer_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_unit_centro_monoerV2.bed'
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
    main(unit_intact_anno_file,unit_monomer_anno_file)