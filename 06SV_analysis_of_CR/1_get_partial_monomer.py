'''
Descripttion: file1:monomer注释文件;file2:gene注释文件;筛选在gene注释末端的monomer序列
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-03-17 23:21:31
'''
import sys
import re
import logging

def merge_ranges(ranges):
    if not ranges:
        return []
    # 将范围按照其起始值进行排序
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    merged_ranges = [sorted_ranges[0]]
    for start, end in sorted_ranges[1:]:
        # 检查当前范围与前一个范围的距离
        if start - merged_ranges[-1][1] <= 10000:
            # 合并范围
            merged_ranges[-1] = [merged_ranges[-1][0], max(end, merged_ranges[-1][1])]
        else:
            # 不需要合并，直接添加到结果列表中
            merged_ranges.append([start, end])
    return merged_ranges


def has_intersection(range_list, rangeA):
    # 分别提取两个范围的起始点和结束点
    start1, end1 = rangeA
    c_has_inter=False
    # print(rangeA,range_list)
    for i in range_list:
        if max(start1, i[0]) <= min(end1, i[1]):
            c_has_inter=True
    return c_has_inter

def merge_ranges_within_distance(range_A, range_list, max_distance=5000):
    start_A,end_A=range_A
    # 遍历范围A的每个范围
    c_has_inter=False
    for i in range(0,len(range_list)):
        # print(range_list[i])
        start_B,end_B=range_list[i]
        if abs(start_A - end_B) <= max_distance or abs(start_B - end_A) <= max_distance:
            range_list[i]=[min(start_A, start_B), max(end_A, end_B)]
            c_has_inter=True
    if not c_has_inter:
        range_list.append(range_A)
    return range_list

def main(monomer_anno_file,gene_anno_file,centro_region_file,monomer_outbed):
    centro_loc={}
    of=open(monomer_outbed,'w')
    with open(centro_region_file,'r') as cf:
        for i in cf:
            if 'start' in i:
                continue
            line=i.strip().split()
            centro_loc[line[0]+'_centro']=[int(line[1]),int(line[2])]
    
    centro_gene_loc={}
    with open(gene_anno_file,'r') as gf:
        for i in gf:
            if 'start' in i:
                continue
            line=i.strip().split()
            gene_len=int(line[2])-int(line[1])
            centro_len=centro_loc[line[0]][1]-centro_loc[line[0]][0]
            if gene_len > 0.3*centro_len:
                continue
            if line[0] in centro_gene_loc:
                cur_gene_loc=[int(line[1]),int(line[2])]
                centro_gene_loc[line[0]]=merge_ranges_within_distance(cur_gene_loc,centro_gene_loc[line[0]])
            else:
                centro_gene_loc[line[0]]=[[int(line[1]),int(line[2])]]
    centro_gene_loc_maxium={}
    for i in centro_gene_loc:
        centro_gene_loc[i]=merge_ranges(centro_gene_loc[i])
        centro_gene_loc_maxium[i]=max([x[1] for x in centro_gene_loc[i]])

    with open(monomer_anno_file,'r') as mf:
        for i in mf:
            line=i.strip().split()
            centro_start=centro_loc[line[0]+'_centro'][0]
            # print(centro_start)
            if line[0]+'_centro' not in centro_gene_loc:
                # print(len([i for i in centro_gene_loc if 'Chr08' in i]))
                continue
            cur_centro_gene_loc=centro_gene_loc[line[0]+'_centro']
            monomer_start=int(line[1])-centro_start
            monomer_end=int(line[2])-centro_start
            # print(cur_centro_gene_loc)
            # print(monomer_start,monomer_end)
            if has_intersection(cur_centro_gene_loc,[monomer_start,monomer_end]):
                outline='\t'.join([line[0],line[1],line[2],line[0]+'_'+line[4]+'_'+line[7]+'_inner',line[3],line[5]])
                of.write(outline+'\n')
            else:
                if monomer_start > centro_gene_loc_maxium[line[0]+'_centro']:
                    outline='\t'.join([line[0],line[1],line[2],line[0]+'_'+line[4]+'_'+line[7]+'_end',line[3],line[5]])
                    of.write(outline+'\n')
    of.close()

if __name__ == "__main__":
    # monomer_anno_file=sys.argv[1]
    monomer_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\monomer_test.txt'
    #这里的monomer test是只提取了Chr08的monomer
    # gene_anno_file=sys.argv[2]
    gene_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_asm_70_CR_Msu7_clean.bed'
    # centro_region_file=sys.argv[3]
    centro_region_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    monomer_outbed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\Chr08_monomer_besideG.bed'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(monomer_anno_file,gene_anno_file,centro_region_file,monomer_outbed)
