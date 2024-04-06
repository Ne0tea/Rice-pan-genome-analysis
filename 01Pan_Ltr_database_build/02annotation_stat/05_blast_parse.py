'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-08 20:14:34
LastEditors: Ne0tea
LastEditTime: 2023-11-13 21:13:39
'''
import sys


def main(blast,sp1,sp2):
    seq1_dic={}
    seq2_dic={}
    out=open(blast+'_parse','w')
    with open(sp1,'r') as sp1f:
        for i in sp1f:
            line=i.strip().split()
            seq1_dic[line[0]]=line[1]
    with open(sp2,'r') as sp2f:
        for i in sp2f:
            line=i.strip().split()
            seq2_dic[line[0]]=line[1]
    with open(blast,'r') as bf:
        blast_pair={}
        for i in bf:
            line=i.strip().split()
            if line[0] == line[1]:
                continue
            identity=float(line[2])
            align_len=int(line[3])
            e_value=float(line[10])
            bitscore=float(line[11])
            if line[0] not in blast_pair:
                blast_pair[line[0]]=[line[1],e_value,align_len,identity,bitscore]
            else:
                ref,c_e,c_align,c_iden,c_bitscore=blast_pair[line[0]]
                if e_value <= c_e:
                    if align_len > 0.98*c_align:
                        if identity > c_iden:
                            if bitscore > 0.98*c_bitscore:
                                blast_pair[line[0]]=[line[1],e_value,align_len,identity,bitscore]
        for i in blast_pair:
            ref=blast_pair[i][0]
            len1=int(seq1_dic[i])
            len2=int(seq2_dic[ref])
            # len_value=max(len1,len2)
            len_value=len1
            ratio=blast_pair[i][2]/len_value
            line_sen=str(i)+'\t'+str(ref)+'\t'+str(blast_pair[i][3])+'\t'+str(blast_pair[i][2])+'\t'\
                +str(len_value)+'\t'+str(ratio)+'\t'+str(blast_pair[i][1])+'\t'+str(blast_pair[i][4])+'\n'
            out.write(line_sen)
    out.close()

if __name__ == "__main__":
    blast_file=sys.argv[1]
    # samtools faidx CX3_consensi.fa.classified > CX3_consensi.fa.classified.fai
    sp1_file=sys.argv[2]
    # samtools faidx CX20_consensi.fa.classified > CX20_consensi.fa.classified.fai
    sp2_file=sys.argv[3]
    # blast_file=r""
    # sp1_file=r""
    # sp2_file=r""
    
    main(blast_file,sp1_file,sp2_file)