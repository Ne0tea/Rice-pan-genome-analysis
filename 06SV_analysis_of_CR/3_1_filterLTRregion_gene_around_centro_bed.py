'''
Descripttion: gene_file:gene注释的bed文件
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-04-06 15:10:23
'''
import sys
import re
import logging
from _00_utilize import has_intersection,list_range_include


def restore_record(ori_list,new_record,centro_dic):
    '''
    restore 20 seq before current into 'ori_list'
    '''
    pattern=re.compile(r'LOC_Os(.*?)g')#LOC_Os10g
    gene_pattern=re.compile(r'ID=(.*?)\.')
    line=new_record.split()
    c_Chr=line[0].split('_')[0]
    c_centro_start=centro_dic[line[0]][0]
    trans_start=str(int(line[3])-c_centro_start)
    trans_end=str(int(line[4])-c_centro_start)
    c_strand=line[6]
    gene_name=gene_pattern.findall(line[8])[0]
    outline=False
    outline='\t'.join([line[0]+'_centro',trans_start,trans_end,line[0]+'_centro|'+gene_name,c_strand,'upstream'])
    if not outline:
        return ori_list
    if len(ori_list) < 20 :
        ori_list.append(outline)
    else:
        ori_list.append(outline)
        ori_list = ori_list[1:]
    return ori_list

def main(gff3_file,centro_file,ltr_region_file,outfile):
    of=open(outfile,'w')
    pattern=re.compile(r'LOC_Os(.*?)g')#LOC_Os10g
    gene_pattern=re.compile(r'ID=(.*?)\.')
    Sygene_pattern=re.compile(r'ID=(.*?);Name')
    CR_record_number={}
    centro_dic={}
    with open(centro_file,'r') as cf:
        for i in cf:
            line=i.strip().split()
            centro_dic[line[0]]=[int(line[1]),int(line[2])]
            CR_record_number[line[0]]=20
    ltr_region_dic={}
    '''
    储存ltr区域信息
    ltr_region_dic={'Chr01_13-65':[[1,10],[14,18]...],
                    'Chr02_13-65':[[1,10],[14,18]...],}
    '''
    with open(ltr_region_file,'r') as ltrf:
        for i in ltrf:
            line=i.strip().split()
            c_start=int(line[1])
            c_end=int(line[2])
            if line[0] in ltr_region_dic:
                ltr_region_dic[line[0]].append([c_start,c_end])
            else:
                ltr_region_dic[line[0]]=[[c_start,c_end]]
            # if line[0]=='Chr12_13-65':
            #     count+=1
    # for i in ltr_region_dic:
    #     print(len(ltr_region_dic[i]))
    # print(count)

    CR_status=False
    ori_record_list=[]
    with open(gff3_file,'r') as mb:
        last_chr=''
        l_start,l_end=0,0
        for i in mb:
            line=i.strip().split()
            if last_chr!=line[0]:
                print(last_chr)
                ori_record_list=[]
                CR_status=False
                last_chr=line[0]
            c_Chr=line[0].split('_')[0]
            c_ltr_region=ltr_region_dic[line[0]]
            c_start=int(line[3])
            c_end=int(line[4])
            ltr_intersect=list_range_include(c_ltr_region,[c_start,c_end])

            c_strand=line[6]
            c_centro_start=centro_dic[line[0]][0]
            trans_start=str(int(line[3])-c_centro_start)
            trans_end=str(int(line[4])-c_centro_start)
            gene_name=gene_pattern.findall(line[8])[0]

            if ltr_intersect:
                continue
            if 'LOC' in line[8]:
                seq_Chr=pattern.findall(line[8])[0]
                if c_Chr!='Chr'+seq_Chr:
                    continue
            else:
                gene_name=Sygene_pattern.findall(line[8])[0]

            if has_intersection([c_start,c_end],centro_dic[line[0]]):
                if not CR_status:
                    for outline in ori_record_list:
                        of.write(outline+'\n')
                CR_status=True
                outline='\t'.join([line[0]+'_centro',trans_start,trans_end,line[0]+'_centro|'+gene_name,c_strand,'centro'])
                of.write(outline+'\n')
            elif min(l_start,l_end) < centro_dic[line[0]][0] and min(c_start,c_end) > centro_dic[line[0]][1]:
                if not CR_status:
                    for outline in ori_record_list:
                        of.write(outline+'\n')
                    CR_status=True
                if 'Chr' not in gene_name:
                    outline='\t'.join([line[0]+'_centro',trans_start,trans_end,line[0]+'_centro|'+gene_name,c_strand,'downstream'])
                    of.write(outline+'\n')
                    CR_record_number[line[0]]-=1
            else:
                if 'Chr' not in gene_name:
                    ori_record_list=restore_record(ori_record_list,i.strip(),centro_dic)
                    if CR_status and CR_record_number[line[0]]>0:
                        outline='\t'.join([line[0]+'_centro',trans_start,trans_end,line[0]+'_centro|'+gene_name,c_strand,'downstream'])
                        of.write(outline+'\n')
                        CR_record_number[line[0]]-=1
                    elif CR_status and CR_record_number[line[0]]<=0:
                        CR_status=False
            l_start,l_end=int(line[3]),int(line[4])
    of.close()

if __name__ == "__main__":
    gff3_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_0318_sorted_filter.gff3'
    # gff3_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_0318_sorted.gff3'
    centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    ltr_region_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_asm_ltr_region.bed'# generate from 3_0_merge_ltr_region.py
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_CR_20.filter.bed'

    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(gff3_file,centro_file,ltr_region_file,outfile)
