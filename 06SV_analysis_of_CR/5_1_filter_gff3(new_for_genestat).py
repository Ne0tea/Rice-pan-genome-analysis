'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-03-21 21:25:05
LastEditors: Ne0tea
LastEditTime: 2024-04-05 12:11:21
'''
import os
import re
import pandas as pd
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
from _00_utilize import has_intersection,list_range_include,intersection_length
# 解析GFF3文件
def parse_gff3(gff_file):
    with open(gff_file) as handle:
        for rec in GFF.parse(handle):
            for feature in rec.features:
                if feature.type == "gene":
                    gene_id = rec.id
                    exons = [exon for mRNA in feature.sub_features for exon in mRNA.sub_features if exon.type == "exon"]
                    yield feature,gene_id, exons

# 过滤基因，如果一个gene条目中包含的exon之间的距离大于2000bp，则过滤掉这个基因
def filter_genes(genes, max_distance=2000):
    for feature,gene_id, exons in genes:
        if len(exons) > 1:
            # 对exon按照起始位置排序
            exons.sort(key=lambda x: x.location.start)
            # 检查相邻exon之间的距离是否大于2000bp
            distances = [exons[i+1].location.start - exons[i].location.end for i in range(len(exons)-1)]
            # print(distances)
            if all(distance <= max_distance for distance in distances):
                yield feature,gene_id, exons
        else:
            yield feature,gene_id, exons

def write_filtered_gff3(filtered_genes, output_gff_file):
    with open(output_gff_file, 'w') as out_handle:
        for gene_feature,gene_id, exons in filtered_genes:
            gene = SeqRecord(seq='', id=gene_id, description="")
            # print(gene_id)
            top_feature=gene_feature
            top_feature.sub_features=exons
            gene.features = [top_feature]
            GFF.write([gene], out_handle)

def read_region_info(centro_file,ltr_region_file):
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
    return centro_dic,ltr_region_dic

def main(gff_file,centro_file,ltr_region_file):
    pattern=re.compile(r'LOC_Os(.*?)g')#LOC_Os10g
    gene_pattern=re.compile(r'ID=(.*?)\.')
    Sygene_pattern=re.compile(r'ID=(.*?);Name')

    tmp_gff3_filtered=gff_file.rstrip('.gff3') + '_tmp0.gff3'
    if not os.path.exists(tmp_gff3_filtered):
        genes = parse_gff3(gff_file)
        # 过滤基因
        filtered_genes = filter_genes(genes)
        # 输出过滤结果
        write_filtered_gff3(filtered_genes, tmp_gff3_filtered)
    tmp_df=pd.read_table(tmp_gff3_filtered,names=['Chr','msu','type','Start','End','prob','Strand','no','ID'])
    tmp_df.sort_values(by=['Chr', 'Start'], inplace=True, ascending=True)
    tmp_df.to_csv(gff_file.rstrip('.gff3') + '_tmp.gff3',sep='\t',index=False,header=False)
    centro_dic,ltr_region_dic=read_region_info(centro_file,ltr_region_file)

    of=open(outfile,'w')
    with open(gff_file.rstrip('.gff3') + '_tmp.gff3','r') as mb:
        last_chr=''
        last_s,last_e=0,0
        for i in mb:
            if i.startswith('#') :
                continue
            line=i.strip().split()
            if line[2] != 'gene':
                continue
            if last_chr!=line[0]:
                dup_gene_name=[]
                print(last_chr)
                last_chr=line[0]
            c_Chr=line[0].split('_')[0]
            ID_rm_CR=line[0].replace('_centro', "")
            c_ltr_region=ltr_region_dic[ID_rm_CR]
            c_centro_start=centro_dic[ID_rm_CR][0]
            c_start=int(float(line[3]))+c_centro_start
            c_end=int(float(line[4]))+c_centro_start
            ltr_intersect=list_range_include(c_ltr_region,[c_start,c_end])

            if ltr_intersect:
                continue
            if 'LOC' in line[8]:
                gene_name=gene_pattern.findall(line[8])[0]
            else:
                gene_name=Sygene_pattern.findall(line[8])[0]
            line=i.strip().split()
            c_Chr=line[0].split('_')[0]
            c_s=c_start
            e_s=c_end
            # if e_s - c_s > 50000:
            #     continue
            ratio=round(intersection_length([last_s,last_e],[c_s,e_s])/(e_s-c_s),4)
            if ratio > 0.9:
                dup_gene_name.append(gene_name)
                continue
            else:
                if c_s <= last_s and e_s >= last_e:
                    dup_gene_name.append(gene_name)
                    continue
            if has_intersection([c_start,c_end],centro_dic[ID_rm_CR]):
                dup_gene_name.append(gene_name)
                of.write(i)
            last_s=c_s
            last_e=e_s
    of.close()

if __name__ == "__main__":
    gff_file = r"E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_asm_70_CR_Msu7_clean.gff3"
    centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    ltr_region_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_asm_ltr_region.bed'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_asm_70_CR_Msu7_clean_filtered2.gff3'
    main(gff_file,centro_file,ltr_region_file)