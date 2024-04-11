'''
Descripttion: 用于gene synteny构建的gff3文件过滤。gene_file:gene注释的bed文件;
标准:   length<150000;
        不允许非同染色体基因注释,允许unscaffold基因注释;
        重叠基因(range dup > 0.9) 不保留;
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-04-02 11:09:51
'''
import sys
import re
import logging

def intersection_length(range1, range2):
    start1, end1 = range1
    start2, end2 = range2
    
    # 计算交集的起始点和结束点
    start = max(start1, start2)
    end = min(end1, end2)
    
    # 如果交集存在，则返回交集长度；否则返回0
    if start <= end:
        return end - start + 1
    else:
        return 0

def main(gff3_file,outfile,centro_file,cross_centro_gene_bed):
    pattern=re.compile(r'LOC_Os(.*?)g')#LOC_Os10g
    cross_centro_gene_of=open(cross_centro_gene_bed,'w')
    centro_dic={}
    with open(centro_file,'r') as cf:
        for i in cf:
            line=i.strip().split()
            centro_dic[line[0]]=[int(line[1]),int(line[2])]
    of=open(outfile,'w')
    with open(gff3_file,'r') as gf:
        last_s=0
        last_e=0
        for i in gf:
            line=i.strip().split()
            c_Chr=line[0].split('_')[0]
            c_s=int(line[3])
            e_s=int(line[4])
            c_centro_region=centro_dic[line[0]]
            if 'LOC' in line[8]:
                seq_Chr=pattern.findall(line[8])[0]
                if c_Chr!='Chr'+seq_Chr:
                    continue
            if c_s < c_centro_region[0] and e_s > c_centro_region[1]:
                cross_centro_gene_of.write(i)
            if e_s - c_s > 150000:
                continue
            ratio=round(intersection_length([last_s,last_e],[c_s,e_s])/(e_s-c_s),4)
            if ratio > 0.9:
                continue
            else:
                of.write(i)
            last_s=c_s
            last_e=e_s
    of.close()
    cross_centro_gene_of.close()


if __name__ == "__main__":
    #generate by sort -V -k 1,1 -k 4,4 all_70_0318.gff3
    gff3_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_0318_sorted.gff3'
    centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_0318_sorted_filter.gff3'
    cross_centro_gene_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\cross_centro_gene.gff3'# gene which over centro region
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(gff3_file,outfile,centro_file,cross_centro_gene_bed)
