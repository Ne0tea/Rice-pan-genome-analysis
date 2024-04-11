'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-12-17 16:08:08
LastEditors: Ne0tea
LastEditTime: 2023-12-20 10:47:51
'''
'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-12-17 16:08:08
LastEditors: Ne0tea
LastEditTime: 2023-12-17 17:09:51
'''
import re
import sys
import os
import matplotlib.pyplot as plt


def main(pan_gff_file_list,pan_centro_bed_file,fai_file,target_seq):
    outfile=open(os.path.join("E:\Bio_analysis\Weedyrice\pan_weedyrice",target_seq+'_len_distribution.txt'),'w')
    outfile_bed=open(os.path.join("E:\Bio_analysis\Weedyrice\pan_weedyrice",target_seq+'_distribution.bed'),'w')
    pattern=re.compile(r'\"Motif:(.*)\"')
    colors=[]
    chr_length={}
    with open(fai_file,'r') as faif:
        for i in faif:
            line=i.strip().split()
            chr_length[line[0]]=int(line[1])
    chr_centro_bed_dic={}
    with open(pan_centro_bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            if line[0].split('_')[0] in chr_centro_bed_dic:
                chr_centro_bed_dic[line[0].split('_')[0]][line[0].split('_')[1]]=[int(line[1]),int(line[2])]
            else:
                chr_centro_bed_dic[line[0].split('_')[0]]={line[0].split('_')[1]:[int(line[1]),int(line[2])]}
    gff_list=[]
    with open(pan_gff_file_list,'r') as gf:
        for i in gf:
            line=i.strip()
            gff_list.append(line)
    point_x=[]#seq in Out centro
    point_y=[]#seq in centro
    for file in gff_list:
        sp=os.path.basename(file).split('.')[0].replace('rmfup_replaced','')
        centro_length={}
        outcentro_length={}
        with open(file,'r') as gffile:
            for i in gffile:
                line=i.strip().split()
                seq_name=re.sub(r'_OS|-I|-LTR|_LTR|_I','',pattern.findall(line[9])[0])
                if seq_name==target_seq:
                    chr=line[0].split('_')[0]
                    sp=line[0].split('_')[1]
                    length=int(line[4])-int(line[3])
                    if (chr=='Chr07' and sp=='CW06') or int(line[3]) > chr_centro_bed_dic[chr][sp][1] \
                        or int(line[4]) < chr_centro_bed_dic[chr][sp][0] :
                        if chr in outcentro_length:
                            outcentro_length[chr]+=length
                        else:
                            outcentro_length[chr]=length
                    else:
                        if chr in centro_length:
                            centro_length[chr]+=length
                        else:
                            centro_length[chr]=length
                    line='\t'.join([line[0],str(line[3]),str(line[4]),str(seq_name),'0',str(line[6])])
                    outfile_bed.write(line+'\n')
        for x in centro_length:
            c_chrl=chr_length[x+'_'+sp]
            c_centrol=chr_centro_bed_dic[x][sp][1]-chr_centro_bed_dic[x][sp][0]
            ratio=(c_chrl-c_centrol)/c_centrol
            line='\t'.join([sp,x,str(centro_length[x]),str(outcentro_length[x]),str(c_centrol),str(c_chrl)])
            outfile.write(line+'\n')
    outfile.close()
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff_rmfup.list"
    # pan centro bed
    # Chr01_13-65 16749781 17724890
    # pan_centro_bed_file=sys.argv[2]
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed"
    # rep_id=sys.argv[3]
    rep_id=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\repeatmask_lib.id'
    fai_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\all_asm_70.fa.fai'

    # target_seq=sys.argv[4]
    target_seq='SZ-22'

    main(pan_gff_file_list,pan_centro_bed_file,fai_file,target_seq)