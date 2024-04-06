'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-12 23:57:34
LastEditors: Ne0tea
LastEditTime: 2023-11-29 10:56:59
'''
import sys


def main(bed_file,num_c=10,len_c=5000):
    seq_dic={}
    with open(bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            length=int(line[2])-int(line[1])
            if 'family' in line[3]:
                if line[3] in seq_dic:
                    seq_dic[line[3]][0]+=1
                    seq_dic[line[3]][1]+=length
                else:
                    seq_dic[line[3]]=[1,length]
    with open(bed_file+'conrtibution.txt','w') as of:
        for i in seq_dic:
            seq=i.strip('\"Motif:').strip(':')
            if seq_dic[i][0] > int(num_c) and seq_dic[i][1] > int(len_c):
                # line=seq+"\t"+str(seq_dic[i][0])+"\t"+str(seq_dic[i][1])
                line=seq
                of.write(line+'\n')
if __name__ == "__main__":
    # bed_file=sys.argv[1]
    bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\NIP_repdeep_rep_sub.bed"
    try:
        num_c=sys.argv[2]
        len_c=sys.argv[3]
    except:
        num_c,len_c=10,5000
    main(bed_file,num_c,len_c)