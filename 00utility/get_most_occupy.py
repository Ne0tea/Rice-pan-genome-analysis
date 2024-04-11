import pandas as pd
import re
import os
import sys
from collections import Counter

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

def main(pan_gff_file_list,pan_centro_bed_file,rep_id,centro_state):
    chr_centro_bed_dic={}
    chr_length_dic={}
    Top_5_sp=pd.DataFrame()
    Class_precentage=pd.DataFrame()
    Family_precentage=pd.DataFrame()
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
    seq_class={}
    with open(rep_id,'r') as repf:
        for i in repf:
            line=i.strip().split('#')
            seq_name=re.sub(r'_OS|-I|-LTR|_LTR|_I','',line[0])
            try:
                cur_seq_class=line[1]
            except IndexError:
                continue
            seq_class[seq_name]=cur_seq_class.split('/')
    pattern=re.compile(r'\"Motif:(.*)\"')
    for file in gff_list:
        # target_motif=[]
        sp_top5_df=pd.DataFrame()
        sp_class_df=pd.DataFrame()
        sp_family_df=pd.DataFrame()
        sp=os.path.basename(file).split('.')[0].replace('rmfup_replaced','')
        motif_chr_dic={}
        if centro_state:
            centrofile=open(file+'centro.region','w')
        else:
            centrofile=open(file+'outcentro.region','w')
        with open(file,'r') as gffile:
            count=0
            for i in gffile:
                if i.startswith('#'):
                    continue
                line=i.strip().split()
                length=int(line[4])-int(line[3])
                seq_name=re.sub(r'_OS|-I|-LTR|_LTR|_I','',pattern.findall(line[9])[0])
                chr_name=line[0].split('_')[0]
                if centro_state:
                    if (chr_name=='Chr07' and sp=='CW06') or int(line[3]) > chr_centro_bed_dic[chr_name][sp][1] \
                        or int(line[4]) < chr_centro_bed_dic[chr_name][sp][0] :
                        continue
                else:
                    if (chr_name=='Chr07' and sp=='CW06') or (int(line[4]) < chr_centro_bed_dic[chr_name][sp][1] and int(line[4]) > chr_centro_bed_dic[chr_name][sp][0]) or \
                        (int(line[3]) > chr_centro_bed_dic[chr_name][sp][0] and int(line[3]) < chr_centro_bed_dic[chr_name][sp][1]) :
                        continue
                centrofile.write(i)
                centro_length=chr_centro_bed_dic[chr_name][sp][1]
                # motif_chr_dic=append_dic(chr_name,seq_name,length,motif_chr_dic,only_length=True)
                motif_chr_dic=append_dic(chr_name,seq_name,length,motif_chr_dic)
                count+=1
        centrofile.close()
        for i in motif_chr_dic:
            chr_dic={}
            num=0
            centro_length=int(chr_centro_bed_dic[i][sp][1])-int(chr_centro_bed_dic[i][sp][0])
            # top5_ratio = [int(x)/centro_length for x in sorted(motif_chr_dic[i].values(),reverse=True)[:5]]
            top5_keys = sorted(motif_chr_dic[i], key=lambda x: motif_chr_dic[i][x],reverse=True)[:5]
            for x in top5_keys:
                num+=1
                chr_dic[i+'_'+str(num)]=x+'__'+str(motif_chr_dic[i][x][0])+'/'+str(motif_chr_dic[i][x][1])+'/'+str(round(motif_chr_dic[i][x][1]/centro_length, ndigits=3))
            keys_line=pd.Series(chr_dic,dtype='str').to_frame()
            sp_top5_df=pd.concat([sp_top5_df,keys_line])
            ###centro region class percentage stat
            class_list=[seq_class[x][0] for x in motif_chr_dic[i] if x in seq_class]
            region_ratio={}
            for x in class_list:
                class_len=sum([motif_chr_dic[i][y][1] for y in motif_chr_dic[i] if y in seq_class and seq_class[y][0]==x])
                region_ratio[i+'_'+x]=class_len/centro_length
                # region_ratio[i+'_'+x]="{:.2%}".format(class_len/centro_length)
            element_count = Counter(class_list)
            total_elements = len(class_list)
            # most_common_five = element_count.most_common(5)
            element_count = Counter(class_list)
            total_elements = len(class_list)
            # element_percentages = {i+'_'+key: "{:.2%}".format(value / total_elements) for key, value in element_count.items()}
            # class_line=pd.Series(element_percentages,dtype='str').to_frame()
            count_kv=sorted(region_ratio.items(), key=lambda x: x[1], reverse=True)[:5]
            top5_reigon={key:"{:.2%}".format(value) for key, value in count_kv}
            class_line=pd.Series(top5_reigon,dtype='str').to_frame()
            sp_class_df=pd.concat([sp_class_df,class_line])
            ###centro region family percentage stat
            family_list=[seq_class[x][1] for x in motif_chr_dic[i] if x in seq_class if len(seq_class[x])>1]
            region_ratio={}
            for x in family_list:
                family_len=sum([motif_chr_dic[i][y][1] for y in motif_chr_dic[i] if y in seq_class and len(seq_class[y])>1 and seq_class[y][1]==x])
                region_ratio[i+'_'+x]=family_len/centro_length
                # region_ratio[i+'_'+x]="{:.2%}".format(family_len/centro_length)
            element_count = Counter(family_list)
            total_elements = len(family_list)
            # most_common_five = element_count.most_common(5)
            element_count = Counter(family_list)
            total_elements = len(family_list)
            # element_percentages = {i+'_'+key: "{:.2%}".format(value / total_elements) for key, value in element_count.items()}
            # family_line=pd.Series(element_percentages,dtype='str').to_frame()
            count_kv=sorted(region_ratio.items(), key=lambda x: x[1], reverse=True)[:5]
            top5_reigon={key:"{:.2%}".format(value) for key, value in count_kv}
            family_line=pd.Series(top5_reigon,dtype='str').to_frame()
            sp_family_df=pd.concat([sp_family_df,family_line])
            # print(family_line)

        sp_top5_df=sp_top5_df.rename(columns={0:sp})
        sp_class_df=sp_class_df.rename(columns={0:sp})
        sp_family_df=sp_family_df.rename(columns={0:sp})
        Top_5_sp=pd.concat([Top_5_sp,sp_top5_df],axis=1)
        Class_precentage=pd.concat([Class_precentage,sp_class_df],axis=1)
        Family_precentage=pd.concat([Family_precentage,sp_family_df],axis=1)
        print(sp,'enrich done!')
    Top_5_sp=Top_5_sp.sort_index()
    Class_precentage=Class_precentage.sort_index()
    Family_precentage=Family_precentage.sort_index()
    if centro_state:
        Top_5_sp.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Centro_region_enrich_seq.csv')
        Class_precentage.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Centro_region_enrich_class.csv')
        Family_precentage.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Centro_region_enrich_family.csv')
    else:
        Top_5_sp.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Centro_out_region_enrich.csv')
        Class_precentage.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Centro_out_region_enrich_class.csv')
        Family_precentage.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Centro_out_region_enrich_family.csv')
    # print(Top_5_sp)
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff_rmfup.list"
    # pan centro bed
    # Chr01_13-65 16749781 17724890
    # pan_centro_bed_file=sys.argv[2]
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\70_pan_centro.bed"
    # rep_id=sys.argv[3]
    rep_id=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\repeatmask_lib.id'
    # centro_state=sys.argv[4]
    centro_state=1
    if int(centro_state) == 0:
        centro_stat=False
    else:
        centro_stat=True
    main(pan_gff_file_list,pan_centro_bed_file,rep_id,centro_state)