'''
Descripttion: gene_file:gene注释的bed文件
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-03-21 10:50:52
'''
import sys
import re
import logging

def has_intersection(rangeA, rangeB):
    # 分别提取两个范围的起始点和结束点
    start1, end1 = rangeA
    start2, end2 = rangeB
    c_has_inter=False
    if max(start1, start2) <= min(end1, end2):
        c_has_inter=True
    return c_has_inter

def main(gff3_file,centro_file,outfile):
    of=open(outfile,'w')
    pattern=re.compile(r'LOC_Os(.*?)g')#LOC_Os10g
    gene_pattern=re.compile(r'ID=(.*?);Name')
    centro_dic={}
    count=0
    with open(centro_file,'r') as cf:
        for i in cf:
            line=i.strip().split()
            centro_dic[line[0]]=[int(line[1]),int(line[2])]
    with open(gff3_file,'r') as mb:
        for i in mb:
            line=i.strip().split()
            c_Chr=line[0].split('_')[0]
            c_start=int(line[3])
            c_end=int(line[4])
            c_strand=line[6]
            c_centro_start=centro_dic[line[0]][0]
            # if c_end - c_start > 300000:
            #     print(gene_name)
            #     continue
            if has_intersection([c_start,c_end],centro_dic[line[0]]):
                trans_start=str(int(line[3])-c_centro_start)
                trans_end=str(int(line[4])-c_centro_start)
                gene_name=gene_pattern.findall(line[8])[0]
                if 'LOC' in line[8]:
                    seq_Chr=pattern.findall(line[8])[0]
                    print(c_Chr,seq_Chr)
                    if c_Chr=='Chr'+seq_Chr:
                        count+=1
                        outline='\t'.join([line[0]+'_centro',trans_start,trans_end,line[0]+'_centro'+gene_name,c_strand])
                        of.write(outline+'\n')
                else:
                    outline='\t'.join([line[0]+'_centro',trans_start,trans_end,line[0]+'_centro'+gene_name,c_strand])
                    of.write(outline+'\n')
    print(count)

if __name__ == "__main__":
    gff3_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_0318.gff3'
    centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_0318.filter.bed'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(gff3_file,centro_file,outfile)
