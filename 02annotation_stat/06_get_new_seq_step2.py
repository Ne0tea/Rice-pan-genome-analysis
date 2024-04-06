'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-12 14:04:14
LastEditors: Ne0tea
LastEditTime: 2023-11-13 19:49:14
'''
import sys
import os

def main(parsed_blast_list_file,repeat_seq_file):
    sp_seq_id=[]
    sp=repeat_seq_file.split('_')[0]
    with open(repeat_seq_file,'r') as listf:
        for i in listf:
            sp_seq_id.append(i.strip())
    blast_list=[]
    with open(parsed_blast_list_file,'r') as listf:
        for i in listf:
            blast_list.append(i.strip())
    seq_link={}
    seq_path={}
    rep_seq=[]
    for fdir in blast_list:
        file=os.path.basename(fdir)
        if file.startswith(sp):
            sp_loc=0
        else:
            sp_loc=1
        norep=True
        if 'repbase' in file:
            repbase=True
        else:
            norep=False
            repbase=False
        outfile=open(file+'_same_seq.id','w')
        sp_list=[file.split('_')[0],file.split('_')[4]]
        sp1=sp
        sp2=list(set(sp_list)-set(sp))[0]
        with open(fdir,'r') as pf:
            for i in pf:
                line=i.strip().split()
                identity=float(line[2])
                alig_len=int(line[3])
                ratio=float(line[5])
                if line[sp_loc].split("#")[0] in sp_seq_id:
                    if ratio>0.6 and identity>60:
                        sentence=line[0]+'\t'+line[1]+'\n'
                        outfile.write(sentence)
                        if not repbase:
                            if sp1+'_'+line[sp_loc] in seq_link:
                                seq_link[sp+'_'+line[sp_loc]].append(sp2+'_'+line[0 if sp_loc ==1 else 1])
                            else:
                                seq_link[sp+'_'+line[sp_loc]]=[sp2+'_'+line[0 if sp_loc ==1 else 1]]
                        else:
                            rep_seq.append(sp1+'_'+line[sp_loc])
        outfile.close()
    merge=list(set(rep_seq) | set(seq_link.keys()))
    # print(len(rep_seq),len(seq_link),len(merge))
    seq_in_repbase=open(sp+'_seq_in_repbase.txt','w')
    seq_local=open(sp+'_seq_local.txt','w')
    seq_in_other_sp=open(sp+'_seq_in_other_sp.txt','w')
    for i in sp_seq_id:
        if sp+'_'+i in [x.split('#')[0] for x in rep_seq]:
            seq_in_repbase.write(i+'\n')
        elif sp+'_'+i in [x.split('#')[0] for x in seq_link]:
            seq_in_other_sp.write(i+'\n')
        else:
            seq_local.write(i+'\n')

if __name__ == "__main__":
    # parsed blast list
    parsed_blast_list_file=sys.argv[1]
    # filtered repeat seq
    species_seq_file=sys.argv[2]
    main(parsed_blast_list_file,species_seq_file)