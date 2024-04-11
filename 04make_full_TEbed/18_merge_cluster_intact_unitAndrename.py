'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-25 20:09:19
LastEditors: Ne0tea
LastEditTime: 2024-02-02 09:25:19
'''
import os
import sys
import re
import copy

def rename_cluster_in_intact(cluster_anno_file,chr_centro_bed_dic,dddidc,unit_name):
    cluster_ord_dic={}
    cur_dic=copy.deepcopy(dddidc)
    for i in list(cur_dic.keys()):
        cur_chr=cur_dic[i][0]
        if unit_name not in cur_dic[i][1]:
            del cur_dic[i]
    with open(cluster_anno_file,'r') as clusterf:
        for i in clusterf:
            line=i.strip().split()
            chr=line[0]
            seq_name=line[3]
            if chr != cur_chr or seq_name != unit_name:
                continue
            # print(line)
            cur_range=str(line[1])+'-'+str(line[2])
            ord=line[5]
            cluster_name=line[6]
            if cluster_name in cluster_ord_dic:
                cluster_ord_dic[cluster_name]+=1
            else:
                cluster_ord_dic[cluster_name]=1
            if cur_range in cur_dic:
                seq_name=unit_name
                # count_ord=cur_dic[cur_range][1].split('_')[6]
                loc=cur_dic[cur_range][1].split('_')[4]
                cur_dic[cur_range][1]=line[0]+'_'+seq_name+'_'+loc+'_'+cluster_name+'_'+'C'+'_'+str(cluster_ord_dic[cluster_name])
            else:
                if int(line[1]) > int(chr_centro_bed_dic[chr][0]) and int(line[2]) < int(chr_centro_bed_dic[chr][1]):
                    loc='Inperi'
                else:
                    loc='Outperi'
                seq_name=unit_name
                # cur_dic[cur_range]=[line[0],seq_name+'_'+loc+'_'+cluster_name+'_'+str(culster_count),ord]
                cur_dic[cur_range]=[line[0],line[0]+'_'+seq_name+'_'+loc+'_'+cluster_name+'_'+'C'+'_'+str(cluster_ord_dic[cluster_name]),ord]
    # print(culster_count)
    return cur_dic

def main(intact_anno_file,cluster_anno_file,pan_centro_bed_file,unit_name):
    # Of=open(Out_file,'w')
    out_dir=os.path.dirname(intact_anno_file)
    merge_dic={}
    pattern=re.compile(r'\"Motif:(.*)\"')
    count=0
    chr_centro_bed_dic={}
    with open(pan_centro_bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            chr_centro_bed_dic[line[0]]=[int(line[1]),int(line[2])]
    intact_process=[]
    perfix='_'.join(os.path.basename(intact_anno_file).split('_')[2:4]).strip('.bed')
    intact_name=None
    with open(intact_anno_file,'r') as intactf:
        for i in intactf:
            line=i.strip().split()
            seq_name=pattern.findall(line[14])[0]
            ord=line[5]
            loc=line[6]
            if unit_name in seq_name:
                if intact_name != line[3]:
                    count+=1
                else:
                    if R_ord == '2':
                        del merge_dic[last_seq_loc]
                if 'LTR' in seq_name:
                    if line[3] in intact_process :
                        R_ord='2'
                    else:
                        R_ord='1'
                        intact_process.append(line[3])
                merge_dic[str(line[1])+'-'+str(line[2])]=[line[0],line[0]+'_'+seq_name+'_'+loc+'_intact_'+str(count)+'_'+R_ord,ord]
                intact_name=line[3]
                last_seq_loc=str(line[1])+'-'+str(line[2])
    Of_I=open(os.path.join(out_dir,'merge_'+perfix+'_'+unit_name+'_I.bed'),'w')
    if unit_name == 'SZ-22':
        merge_I_dic=rename_cluster_in_intact(cluster_anno_file,chr_centro_bed_dic,merge_dic,'SZ-22_I')
    else:
        merge_I_dic=copy.deepcopy(merge_dic)
        for i in list(merge_I_dic.keys()):
            if unit_name+'_I' not in merge_I_dic[i][1]:
                del merge_I_dic[i]
                continue
    for i in merge_I_dic:
        line=merge_I_dic[i][0]+'\t'+'\t'.join(i.split('-'))+'\t'+merge_I_dic[i][1]+'\t'+'0'+'\t'+merge_I_dic[i][2]
        Of_I.write(line+'\n')
    Of_I.close()

    Of_LTR=open(os.path.join(out_dir,'merge_'+perfix+'_'+unit_name+'_LTR.bed'),'w',encoding='UTF-8')
    if unit_name == 'SZ-22':
        merge_LTR_dic=rename_cluster_in_intact(cluster_anno_file,chr_centro_bed_dic,merge_dic,'SZ-22_LTR')
    else:
        merge_LTR_dic=copy.deepcopy(merge_dic)
        for i in list(merge_LTR_dic.keys()):
            if unit_name+'_LTR' not in merge_LTR_dic[i][1]:
                del merge_LTR_dic[i]
                continue
    for i in merge_LTR_dic:
        line=merge_LTR_dic[i][0]+'\t'+'\t'.join(i.split('-'))+'\t'+merge_LTR_dic[i][1]+'\t'+'0'+'\t'+merge_LTR_dic[i][2]
        Of_LTR.write(line+'\n')
    Of_LTR.close()

if __name__ == "__main__":
    #python 18_merge_cluster_intact_unitAndrename.py intact_anno_file cluster_anno_file pan_centro_bed_file SZ-22
    # intact_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\intact_unit_Chr08_NJ11.bed'
    intact_anno_file=sys.argv[1]
    # cluster_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22_cluster_distribution.bed'
    cluster_anno_file=sys.argv[2]
    # pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed"
    pan_centro_bed_file=sys.argv[3]
    # unit_name='SZ-22'
    unit_name=sys.argv[4]
    main(intact_anno_file,cluster_anno_file,pan_centro_bed_file,unit_name)