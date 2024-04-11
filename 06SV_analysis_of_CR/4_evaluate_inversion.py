'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-03-28 16:46:43
LastEditors: Ne0tea
LastEditTime: 2024-03-30 13:13:38
'''
import logging
from Bio.Align import PairwiseAligner
from Bio import pairwise2
import numpy as np
import pandas as pd
from Bio.pairwise2 import format_alignment

def DP_align_backup(seq1,seq2):
    if len(seq1) <= len(seq2):
        seq1=seq1[::-1]
    else:
        seq2=seq2[::-1]
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -2
    alignments = aligner.align(seq1,seq2)
    alignments_score=[alignment.score for alignment in sorted(alignments)]
    return max(alignments_score)

def DP_align(seq1,seq2,score_only=False):
    if len(seq1) <= len(seq2):
        seq1=seq1[::-1]
        full_match_score=len(seq1)*2
    else:
        seq2=seq2[::-1]
        full_match_score=len(seq2)*2
    if score_only:
        max_score=pairwise2.align.globalms(seq1, seq2,match=2, mismatch=-1, open=-3,extend=-2,score_only=True)
        return max_score/full_match_score
    else:
        align_result=pairwise2.align.globalms(seq1, seq2,match=2, mismatch=-1, open=-3,extend=-2)
        max_score=max([x.score for x in align_result])
        return align_result,max_score/full_match_score

def merge_adjacent_duplicates(lst):
    merged_list = []
    i = 0
    while i < len(lst):
        current_element = lst[i]
        merged_list.append(current_element)
        while i + 1 < len(lst) and lst[i + 1] == current_element:
            i += 1
        i += 1
    return merged_list

def main(monomer_bed_file,phy_ord_file,score_only):
    convert_keydic='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    phy_ord=[]
    with open(phy_ord_file,'r') as phyf:
        for i in phyf:
            line=i.strip().split()
            phy_ord.append(line[0])
    with open(monomer_bed_file,'r') as mfile:
        monomer_record_dic={}
        for i in mfile:
            line=i.strip().split()
            cur_chr,cur_sp=line[0].split('_')[:]
            if cur_chr in monomer_record_dic:
                if cur_sp in monomer_record_dic[cur_chr]:
                    monomer_record_dic[cur_chr][cur_sp].append(line[5])
                else:
                    monomer_record_dic[cur_chr][cur_sp]=[line[5]]
            else:
                monomer_record_dic[cur_chr]={cur_sp:[line[5]]}
    unique_monomer_record_dic={}
    sp_list={}
    for cur_chr in monomer_record_dic:
        sp_list[cur_chr]=[]
        for cur_sp in monomer_record_dic[cur_chr]:
            sp_list[cur_chr].append(cur_sp)
            cur_chr_sp_monomer_list=monomer_record_dic[cur_chr][cur_sp]
            unique_cur_chr_sp_monomer_list=merge_adjacent_duplicates(cur_chr_sp_monomer_list)
            if cur_chr in unique_monomer_record_dic:
                    unique_monomer_record_dic[cur_chr].append(unique_cur_chr_sp_monomer_list)
            else:
                unique_monomer_record_dic[cur_chr]=[unique_cur_chr_sp_monomer_list]
            # if cur_sp=='CW15':
            #     print(unique_cur_chr_sp_monomer_list)
    del monomer_record_dic
    for chr in unique_monomer_record_dic:
        cur_order=sp_list[chr]
        cur_chr_align_score_matrix=np.zeros((len(cur_order), len(cur_order)))
        for sp1 in list(range(0,len(unique_monomer_record_dic[chr])-1)):
            tmp_seq1=unique_monomer_record_dic[chr][sp1]
            for sp2 in list(range(sp1+1,len(unique_monomer_record_dic[chr]))):
                tmp_seq2=unique_monomer_record_dic[chr][sp2]
                transvert_ori=list(set(tmp_seq1+tmp_seq2))
                transvert_dic={transvert_ori[x]:convert_keydic[x] for x in range(len(transvert_ori))}
                sp1_seq=''.join([transvert_dic[x] for x in tmp_seq1])
                sp2_seq=''.join([transvert_dic[x] for x in tmp_seq2])
                sp1_seq=sp1_seq[:50]+sp1_seq[-50:]
                sp2_seq=sp2_seq[:50]+sp2_seq[-50:]
                if score_only:
                    cur_score=DP_align(sp1_seq,sp2_seq,score_only)
                else:
                    align_result,cur_score=DP_align(sp1_seq,sp2_seq)
                    if cur_order[sp1]=='13-65' and cur_order[sp2]=='HuaZhan':
                        print(tmp_seq1,tmp_seq2)
                        print(format_alignment(*align_result[0]))
                        print(cur_score)
                cur_chr_align_score_matrix[sp1,sp2]=cur_score
        lower_triangle_indices = np.tril_indices_from(cur_chr_align_score_matrix, -1)
        cur_chr_align_score_matrix[lower_triangle_indices] = cur_chr_align_score_matrix.T[lower_triangle_indices]

        df = pd.DataFrame(cur_chr_align_score_matrix)
        df.index = cur_order  # 行名
        df.columns = cur_order  # 列名
        # 使用reindex()方法排序
        df_sorted = df.reindex(index=phy_ord, columns=phy_ord)
        df_sorted.to_csv('matrix.csv')
    print(cur_chr_align_score_matrix)
if __name__ == "__main__":
    monomer_bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\monomer_test.txt'
    phy_ord_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_phylogenic.txt'
    score_only=True
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(monomer_bed_file,phy_ord_file,score_only)