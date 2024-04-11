'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-12-02 22:00:18
LastEditors: Ne0tea
LastEditTime: 2023-12-02 23:35:15
'''
import sys
hmm_list_file=sys.argv[1]
hmm_list=[]

pf_prot_dic={
    'PF00026.27':'AP','PF13650.10':'AP','PF08284.15':'AP',
    'PF13975.10':'AP','PF00077.24':'AP','PF09668.14':'AP','PF03732.21':'GAG',
    'PF00665.30':'IN','PF13683.10':'IN','PF17921.5':'IN',
    'PF02022.23':'IN','PF09337.14':'IN','PF00552.25':'IN','PF13456.10':'RNase',
    'PF17917.5':'RNase','PF17919.5':'RNase','PF00078.31':'RT'
}
with open(hmm_list_file,'r') as f_list:
    for i in f_list:
        line=i.strip()
        hmm_list.append(line)
seq_dic={}
for i in hmm_list:
    with open(i,'r') as hmmf:
        for x in hmmf:
            if x.startswith('#'):
                continue
            line=x.strip().split()
            seq_name=line[0]
            pfam_acc=line[2]
            pfam_num=line[3]
            e_value=float(line[4])
            class_pfam=pf_prot_dic[pfam_num]
            if e_value < 1:
                if seq_name in seq_dic:
                    seq_dic[seq_name].append(class_pfam)
                else:
                    seq_dic[seq_name]=[class_pfam]
for i in seq_dic:
    if len(set(seq_dic[i])) > 3:
        print(i,'\t'.join(sorted(list(set(seq_dic[i])))),sep='\t')