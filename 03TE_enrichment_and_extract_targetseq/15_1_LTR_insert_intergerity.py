'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-02-12 21:51:12
'''
import sys
import logging
import re
import os
import copy

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

def main(LTR_range_file,fasta_file,unit_name):
    out_file=open(unit_name+'_fa300.bed','w')

    with open(LTR_range_file,'r') as lf:
        tmp_dic={}
        for  i in lf:   
            line=i.strip().split()
            cur_start=int(line[1])
            cur_end=int(line[2])
            unit_ID=line[3].rsplit('_',1)[0]
            if unit_ID in tmp_dic:
                new_start=min(tmp_dic[unit_ID][1],cur_start)
                new_end=max(tmp_dic[unit_ID][2],cur_end)
                tmp_dic[unit_ID][1]=new_start
                tmp_dic[unit_ID][2]=new_end
            else:
                tmp_dic[unit_ID]=[line[0],cur_start,cur_end]
    
    with open('merge_'+unit_name+'.bed','w') as tmpf:
        for i in tmp_dic:
            line=tmp_dic[i][0]+'\t'+str(tmp_dic[i][1])+'\t'+str(tmp_dic[i][2])+'\t'+i
            tmpf.write(line+'\n')
    
    with open('tmp.bed','r') as lf:
        for i in lf:
            line=i.strip().split()
            if unit_name in line[3]:
                # print(1111)
                cur_start=copy.deepcopy(int(line[1]))
                cur_end=copy.deepcopy(int(line[2]))
                cur_name=copy.deepcopy(line[3])
                line[1]=cur_start-300
                line[2]=cur_start
                line[3]=cur_name+'_fw'
                out_line='\t'.join([str(x) for x in line])
                out_file.write(out_line+'\n')
                line[1]=cur_end
                line[2]=cur_end+300
                line[3]=cur_name+'_aw'
                out_line='\t'.join([str(x) for x in line])
                out_file.write(out_line+'\n')
    out_file.close()
    cmd1='bedtools getfasta -name -fi '+fasta_file+' -fo SZ-22_extract_fa300.fasta -bed '+unit_name+'_fa300.bed'
    os.system(cmd1)
    print('Get forward and afterward sequence !')
    cmd2='blastn -task blastn -db Monomer_db -num_threads 10 -query SZ-22_extract_fa300.fasta -outfmt 6 -out SZ-22_extract_fa300.blast_out '
    os.system(cmd2)
    print('Blast finished !')

if __name__ == "__main__":
    # monomer_file=sys.argv[1]#SZ-22_unit_centro.bed 
    # monomer_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\monomer_anno.txt'
    LTR_range_file=sys.argv[1]
    # LTR_range_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\merge_full_LTR.bed'
    fasta_file=sys.argv[2]
    unit_name=sys.argv[3]
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(LTR_range_file,fasta_file,unit_name)
