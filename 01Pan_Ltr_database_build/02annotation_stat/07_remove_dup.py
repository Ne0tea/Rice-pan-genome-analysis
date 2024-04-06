'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-13 21:51:17
LastEditors: Ne0tea
LastEditTime: 2023-11-15 14:03:25
'''
import sys
import re
import os

def main(id_file_list,blast_file_list):
    start_list=[]
    end_list=[]
    cluster_dic={}
    id_trans_dic={}
    redunt_list=[]
    redunt_seq_dic={}
    redunt_id_file=open('redunt_id_file.txt','w')
    with open(id_file_list,'r') as idf:
        for i in idf:
            redunt_list.append(i.strip())
    for file in redunt_list:
        sp=os.path.basename(file)
        sp_name=sp.split('_')[0]
        redunt_seq_dic[sp_name]=[]
        with open(file,'r') as redf:
            for x in redf:
                redunt_id_file.write(sp_name+'_'+x)
                redunt_seq_dic[sp_name].append(sp_name+'_'+x.strip())
    redunt_id_file.close()
    # print(redunt_seq_dic)
    blast_list=[]
    with open(blast_file_list,'r') as blf:
        for i in blf:
            blast_list.append(i.strip())
    for file in blast_list:
        sp=os.path.basename(file)
        sp1_name=sp.split('_')[0]
        sp2_name=sp.split('_')[4]
        with open(file,'r') as spblastfile:
            for i in spblastfile:
                line=re.split('\t|#',i)
                # line=i.strip().split().split('#')
                id_trans_dic[sp1_name+'_'+line[0]]=sp1_name+'_'+line[0]+'#'+line[1]
                id_trans_dic[sp2_name+'_'+line[2]]=sp2_name+'_'+line[0]+'#'+line[3]
    # print(id_trans_dic)
    sp_pairwise_blast=open('sp_pairwise_blast.blastout','w')
    for file in blast_list:
        sp=os.path.basename(file)
        sp1_name=sp.split('_')[0]
        sp2_name=sp.split('_')[4]
        with open(file,'r') as spblastfile:
            for i in spblastfile:
                line=re.split('\t|#',i)
                # line=i.strip().split().split('#')
                # print(sp_name+'_'+line[0],sp1_name)
                if sp1_name+'_'+line[0] in redunt_seq_dic[sp1_name] and \
                    sp2_name+'_'+line[2] in redunt_seq_dic[sp2_name]:
                    sp_pairwise_blast.write(sp1_name+'_'+line[0]+"\t"+sp2_name+'_'+line[2]+'\n')
                    start_list.append(sp1_name+'_'+line[0])
                    end_list.append(sp2_name+'_'+line[2])
    sp_pairwise_blast.close()
    # print(start_list)
    left_search=list(set(start_list))
    # print(len(left_search))
    for cluster_num in range(0,len(start_list)):
        if start_list[cluster_num] not in left_search:
            continue
        c_start_list=list(set(start_list))
        c_end_list=list(set(end_list))
        cluster_dic['cluster'+str(cluster_num)]=[start_list[cluster_num],end_list[cluster_num]]
        # left_search.remove(start_list[cluster_num])

        while len(set(cluster_dic['cluster'+str(cluster_num)]) & set(c_start_list))!=0 :
                    index_s=[ a for x in set(cluster_dic['cluster'+str(cluster_num)]).intersection(c_start_list) for a, b in list(enumerate(start_list)) if b == x]
                    index_e=[a for x in set(cluster_dic['cluster'+str(cluster_num)]).intersection(c_end_list) for a, b in list(enumerate(end_list)) if b == x ]
                    seq_s=[start_list[i] for i in index_e]
                    # print(seq_s)
                    seq_e=[end_list[i] for i in index_s]
                    # print(seq_e)
                    cluster_dic['cluster'+str(cluster_num)].extend(list(set([x for x in seq_s if x not in cluster_dic['cluster'+str(cluster_num)] ])))
                    cluster_dic['cluster'+str(cluster_num)].extend(list(set([x for x in seq_e if x not in cluster_dic['cluster'+str(cluster_num)] ])))
                    # print(len(c_start_list))
                    for x in set(seq_s):
                        c_start_list=list(filter(lambda y:y!=x, c_start_list))
                        # print(len(c_start_list))
                        if x in left_search:
                            left_search.remove(x)
    os.makedirs('Cluster_seq_and_dominate',exist_ok=True)
    os.chdir('./Cluster_seq_and_dominate')
    for i in cluster_dic:
        with open(i+'.id','w') as clusterf:
            for x in cluster_dic[i]:
                clusterf.write(x+'\n')
        with open(i+'true.id','w') as truef:
            true_id=[id_trans_dic[x] for x in cluster_dic[i]]
            for x in true_id:
                truef.write(x.split('#')[0]+'\t'+x+'\n')
        os.system('get_fasta_from_id.pl -s ../all_denovo_lib.fa -f '+i+'.id')
        print('Get '+i+' fasta!')
        os.system('cd-hit -i '+i+'.id.fasta'+' -o '+i+'_cluster'+' -c 0.7 -aS 0.8 -d 0')
        print('Get '+i+'representative fasta!')
    os.chdir('../')
        # print(i,len(set(cluster_dic[i])),cluster_dic[i],sep='\t')
if __name__ == "__main__":
    # redunant id file list
    redunant_id_file_list=sys.argv[1]
    # sp pairwise blast file list
    sp_blast_file_list=sys.argv[2]
    main(redunant_id_file_list,sp_blast_file_list)