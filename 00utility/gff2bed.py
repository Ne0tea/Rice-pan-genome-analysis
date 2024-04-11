'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-26 16:26:50
LastEditors: Ne0tea
LastEditTime: 2023-11-26 17:42:55
'''
import pandas as pd
import re
import pandas as pd
import os
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def append_dic(chr,name,length,dic,only_length=False):
    if chr in dic:
        if name in dic[chr]:
            if only_length:
                dic[chr][name]=dic[chr][name]+length
            else:
                dic[chr][name]=[dic[chr][name][0]+1,dic[chr][name][1]+length]
        else:
            if only_length:
                dic[chr][name]=length
            else:
                dic[chr][name]=[1,length]
    else:
        if only_length:
            dic[chr]={name:length}
        else:
            dic[chr]={name:[1,length]}
    return dic

def main(pan_gff_file_list,pan_centro_bed_file,sp_phy_file,repeat_id,centro_stat):
    chr_centro_bed_dic={}
    chr_length_dic={}
    outfile=open(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_repeat.bed','w')
    uniq_seq=[]
    sp_phy={}
    id_dic={}
    with open(sp_phy_file,'r') as phyf:
        for i in phyf:
            line=i.strip().split()
            sp_phy[line[0]]=line[1]
    # print(sp_phy)
    with open(repeat_id,'r') as repeatf:
        for i in repeatf:
            line=i.strip().split('#')
            if len(line)==1:
                continue
            id_dic[line[0]]=line[1]
    with open(pan_centro_bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            if line[0].split('_')[0] in chr_centro_bed_dic:
                chr_centro_bed_dic[line[0].split('_')[0]][line[0].split('_')[1]]=[int(line[1]),int(line[2])]
            else:
                chr_centro_bed_dic[line[0].split('_')[0]]={line[0].split('_')[1]:[int(line[1]),int(line[2])]}
            if line[0].split('_')[0] in chr_length_dic:
                chr_length_dic[line[0].split('_')[0]][line[0].split('_')[1]]=int(line[3])
            else:
                chr_length_dic[line[0].split('_')[0]]={line[0].split('_')[1]:int(line[3])}
    gff_list=[]
    with open(pan_gff_file_list,'r') as gf:
        for i in gf:
            line=i.strip()
            gff_list.append(line)
    pattern=re.compile(r'\"Motif:(.*)\"')
    for file in gff_list:
        # target_motif=[]
        sp=os.path.basename(file).split('.')[0]
        motif_chr_dic={}
        with open(file,'r') as gffile:
            count=0
            for i in gffile:
                if i.startswith('#'):
                    continue
                line=i.strip().split()
                # length=int(line[4])-int(line[3])
                seq_name=pattern.findall(line[9])[0]
                chr_name=line[0].split('_')[0]
                
                if centro_stat:
                    if (chr_name=='Chr07' and sp=='CW06') or int(line[3]) > chr_centro_bed_dic[chr_name][sp][1] \
                        or int(line[4]) < chr_centro_bed_dic[chr_name][sp][0] :
                        continue
                else:
                    # print(1)
                    if (chr_name=='Chr07' and sp=='CW06') or (int(line[4]) < chr_centro_bed_dic[chr_name][sp][1] and int(line[4]) > chr_centro_bed_dic[chr_name][sp][0]) or \
                        (int(line[3]) > chr_centro_bed_dic[chr_name][sp][0] and int(line[3]) < chr_centro_bed_dic[chr_name][sp][1]) :
                        continue
                centro_length=chr_centro_bed_dic[chr_name][sp][1]
                cur_seq_name=pattern.findall(line[9])[0].split('_')[0]
                if cur_seq_name.endswith('-LTR') :
                    cur_seq_name=cur_seq_name.replace('-LTR','')
                elif cur_seq_name.endswith('-I'):
                    cur_seq_name=cur_seq_name.replace('-I','')
                elif cur_seq_name.endswith('LTR'):
                    cur_seq_name=cur_seq_name.replace('LTR','')
                elif cur_seq_name.endswith('LTRA'):
                    cur_seq_name=cur_seq_name.replace('LTRA','')
                uniq_seq.append(cur_seq_name)
                # for i in set(uniq_seq):
                #     if 'LTR' in i:
                #         print(i)
                # print(line[3]+'\t'+line[4]+'\t'+line[6])
                if seq_name in id_dic:
                    outline=line[0]+'\t'+line[3]+'\t'+line[4]+'\t'+line[6]+'\t'+cur_seq_name+'\t'+id_dic[seq_name]+'\t'+sp+'\t'+sp_phy[sp]+'\n'
                else:
                    outline=line[0]+'\t'+line[3]+'\t'+line[4]+'\t'+line[6]+'\t'+cur_seq_name+'\t'+'Unknow'+'\t'+sp+'\t'+sp_phy[sp]+'\n'
                outfile.write(outline)
                count+=1
    outfile.close()
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff.list"
    # pan centro bed
    # Chr01_13-65 16749781 17724890
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\70_pan_centro.bed"
    sp_phy_file=r'pan_weedyrice/70_phylogenic.txt'
    repeat_id=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\repeatmask_lib.id'
    centro_state=1
    if int(centro_state) == 0:
        centro_stat=False
    else:
        centro_stat=True
    main(pan_gff_file_list,pan_centro_bed_file,sp_phy_file,repeat_id,centro_stat)