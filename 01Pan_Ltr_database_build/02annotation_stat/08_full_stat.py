'''
Descripttion: 根据seq name的特点,区分internal和LTR,将以LTR开始并以重复结束的unit视为intact,反之为unintact
Out:过滤所有染色体上freq<5&&len>7500的类型
Author: Ne0tea
version: 
V1: 将间隔5个其余LTR类型之内的LTR都视为一个完整的LTR
V2: 首尾相连结构(LTR+I+LTR)被视为一个完整的插入结构
V3: 1.限定solo intact 结构元件之间的距离(0.1 concensus seq length)
    2.根据ful元件中插入类型(短串联,卫星序列SAT2),拆分成sat sd hybird类型
V4: 新增fol参数:拆分freq和length的统计
V5: Bug fix:
            1.merge之后的gff文件条目之间的距离变远,导致intact LTR丢失
            2.SAT序列插入LTR类型针对分析,根据重新注释monomer信息分析插入情况
Date: 2023-11-16 19:46:32
LastEditors: Ne0tea
LastEditTime: 2024-01-22 22:18:20
'''
import sys
import re
import pandas as pd
import os
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def criterion_freq(series):
    state=False
    global fol
    # print(series)
    if series.name == 'unit_number':
        return True
    if fol=='freq':
        for i in series:
            if int(i) >= 5 :
                state=True
                return state
    elif fol=='length':
        for i in series:
            if int(i) >= 7500 :
                state=True
                return state
    return state

def criterion_freq2(series):
    state=False
    global fol
    if series.name == 'unit_number':
        return True
    if fol=='freq':
        for i in series:
            if int(i) >= 2 :
                state=True
                return state
    elif fol=='length':
        for i in series:
            if int(i) >= 3000 :
                state=True
                return state
    return state

def append_dic(chr,name,length,dic):
    if chr in dic:
        if name in dic[chr]:
            dic[chr][name]=[dic[chr][name][0]+1,dic[chr][name][1]+length]
        else:
            dic[chr][name]=[1,length]
    else:
        dic[chr]={name:[1,length]}
    return dic

def main(pan_gff_file_list,pan_centro_bed_file,centro_stat,inter_v_ltr:str,ltr_length_file,fol):
    outfile=open('Ful_element.gff3','w')
    ivl_dic={}
    final_df=pd.DataFrame()
    intact_df=pd.DataFrame()
    with open(inter_v_ltr,'r') as ivf:
        for i in ivf:
            line=i.strip().split()
            ivl_dic[line[0]]=line[1]
            ivl_dic[line[1]]=line[0]
    seq_len={}
    with open(ltr_length_file,'r') as lenf:
        for i in lenf:
            line=i.strip().split()
            seq_len[line[0]]=int(line[1])
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
    # print(chr_centro_bed_dic)
    gff_list=[]
    with open(pan_gff_file_list,'r') as gf:
        for i in gf:
            line=i.strip()
            gff_list.append(line)
    pattern=re.compile(r'\"Motif:(.*)\"')
    ful_top5_df=pd.DataFrame()
    if centro_stat:
        ful_unit_bed=open(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\ful_unit_centro.bed','w')
    else:
        ful_unit_bed=open(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\ful_unit_OUT_centro.bed','w')
    # V5 remove
    # if centro_stat:
    #     insert_unit_bed=open(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_insert_unit_centro.bed','w')
    # else:
    #     insert_unit_bed=open(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\ful_insert_unit_OUT_centro.bed','w')
    for file in gff_list:
        target_motif=[]
        sp=os.path.basename(file).split('.')[0].replace('rmfup_replaced','')
        motif_chr_dic={}#{seq_name:[chr_name,start,end,gap_ratio]}
        intact_motif_chr_dic={}#{seq_name:[chr_name,start,end,gap_ratio]}
        with open(file,'r') as gffile:
            count=0
            last_loc=0
            for i in gffile:
                if i.startswith('#'):
                    continue
                line=i.strip().split()
                '''
                length=int(line[4])-int(line[3])
                '''
                start=int(line[3])
                end=int(line[4])
                seq_name=pattern.findall(line[9])[0]
                chr_name=line[0].split('_')[0]
                last_gap=start-last_loc
                ord=line[6]
                if seq_name not in seq_len:
                    seq_concus_len=end-start
                else:
                    seq_concus_len=seq_len[seq_name]
                if centro_stat:
                    if (chr_name=='Chr07' and sp=='CW06') or int(line[3]) > chr_centro_bed_dic[chr_name][sp][1] \
                        or int(line[4]) < chr_centro_bed_dic[chr_name][sp][0] :
                        continue
                else:
                    # print(1)
                    if (chr_name=='Chr07' and sp=='CW06') or (int(line[4]) < chr_centro_bed_dic[chr_name][sp][1] and int(line[4]) > chr_centro_bed_dic[chr_name][sp][0]) or \
                        (int(line[3]) > chr_centro_bed_dic[chr_name][sp][0] and int(line[3]) < chr_centro_bed_dic[chr_name][sp][1]) :
                        continue
                gap_ratio=round(last_gap/seq_concus_len,3)
                intact_ratio=round((end-start)/seq_concus_len,3)
                target_motif.append({seq_name:[chr_name,start,end,gap_ratio,ord,intact_ratio]})
                count+=1
                last_loc=end
        print(sp,'running')
        index=0
        unit_numb={}
        insert_numb={}
        while index < len(target_motif):
            ful_state=True
            seq=list(target_motif[index])[0]
            chr=target_motif[index][seq][0]
            ord=target_motif[index][seq][4]
            if chr in unit_numb:
                unit_numb[chr]+=1
            else:
                unit_numb[chr]=1
            length=target_motif[index][seq][2]-target_motif[index][seq][1]
            '''
            length=target_motif[index][seq][1]
            '''
            if seq in ivl_dic:
                acculate_dic={}
                c_index=index
                if 'LTR' in seq:
                    ltr_name=seq
                    inter_name=ivl_dic[seq]
                else:
                    ful_state=False
                    inter_name=seq
                    ltr_name=ivl_dic[seq]
                if c_index > 0:
                    last_name=list(target_motif[c_index-1])[0]
                    cur_gapR=target_motif[c_index][seq][3]
                if  c_index < len(target_motif)-1:
                    next_name=list(target_motif[c_index+1])[0]
                    next_chr=target_motif[c_index+1][next_name][0]
                    next_gapR=target_motif[c_index+1][next_name][3]
                if ful_state:
                    # print(c_index)
                    if c_index <= len(target_motif)-3 and \
                        (next_name in [inter_name,'SAT-2_OS'] or '(' in next_name) and \
                        chr==next_chr:
                        last_intact=target_motif[c_index][seq][5]
                        cal_index=c_index+1
                        cal_name=list(target_motif[cal_index])[0]
                        cal_chr=target_motif[cal_index][cal_name][0]
                        cal_gapR=target_motif[cal_index][cal_name][3]
                        cal_intact=target_motif[cal_index][cal_name][5]
                        sata=0
                        internal=False
                        unit_list=[seq]
                        cur_ltr_ord='LTR_1'
                        acculate_dic['LTR_1']=last_intact
                        while (cal_name in [ltr_name,inter_name,'SAT-2_OS'] or '(' in cal_name) and\
                                cal_gapR < 0.1 and max(acculate_dic.values()) < 1.2 and\
                                cal_index < len(target_motif)-1 and cal_chr==chr:
                            # if sata>0:
                            #     sata=0
                            # if cal_name=='SAT-2_OS' or '(' in cal_name :
                            #     if last_intact < 0.8:
                            #         sata+=1
                            #         if cal_name=='SAT-2_OS':
                            #             sata_insert=True
                            #         elif '(' in cal_name:
                            #             sd_insert=True
                            #     else:
                            #         break
                            # if cal_name==ltr_name:
                            #     ltr_num+=1
                            if cal_name==inter_name:
                                if 'Inner' in acculate_dic:
                                    acculate_dic['Inner']+=cal_intact
                                else:
                                    acculate_dic['Inner']=cal_intact
                                internal=True
                                cur_ltr_ord='LTR_2'
                                # sata-=1
                            if cal_name==ltr_name:
                                if cur_ltr_ord=='LTR_2' and 'LTR_2' not in acculate_dic:
                                    acculate_dic[cur_ltr_ord]=cal_intact
                                else:
                                    acculate_dic[cur_ltr_ord]+=cal_intact
                            unit_list.append(cal_name)
                            cal_index+=1
                            cal_name=list(target_motif[cal_index])[0]
                            cal_chr=target_motif[cal_index][cal_name][0]
                            cal_gapR=target_motif[cal_index][cal_name][3]
                            cal_intact=target_motif[cal_index][cal_name][5]
                        if max(acculate_dic.values()) >= 1.2 :
                            unit_list.pop()
                            if len(unit_list) < 1:
                                unit_name='nful_'+ltr_name
                                unit_length=sum([int(x[list(x)[0]][2])-int(x[list(x)[0]][1]) for x in target_motif[c_index:c_index+1]])
                                motif_chr_dic=append_dic(chr,unit_name,unit_length,motif_chr_dic)
                                index+=1
                                continue
                        while unit_list[-1]=='SAT-2_OS':
                            unit_list.pop()
                        # elif unit_list[-1] not in [inter_name,'SAT-2_OS'] or '(' in cal_name:
                        #     unit_list.pop()
                        # print(unit_list)
                        if internal and unit_list[-1]==ltr_name:
                            # V5 remove
                            # unit_name='ful_intact_'+ltr_name
                            if chr in insert_numb:
                                insert_numb[chr]+=1
                            else:
                                insert_numb[chr]=1
                            unit_start=target_motif[c_index][seq][1]
                            unit_end=target_motif[c_index+len(unit_list)-1][seq][2]
                            unit_length=sum([int(x[list(x)[0]][2])-int(x[list(x)[0]][1]) for x in target_motif[c_index:c_index+len(unit_list)]])
                            # V5 remove
                            # intact_motif_chr_dic=append_dic(chr,unit_name,unit_length,intact_motif_chr_dic)
                            # if not unit_name.startswith('ful_intact_'):
                            #     line='\t'.join([chr+'_'+sp,str(unit_start),str(unit_end),chr+'_'+sp+'_'+unit_name+str(insert_numb[chr]),'0',ord,'centro' if centro_stat else 'Outcentro'])
                            #     insert_unit_bed.write(line+'\n')
                            unit_name='ful_'+ltr_name
                            motif_chr_dic=append_dic(chr,unit_name,unit_length,motif_chr_dic)
                            line='\t'.join([chr+'_'+sp,str(unit_start),str(unit_end),chr+'_'+sp+'_'+ltr_name+str(insert_numb[chr]),'0',ord,'centro' if centro_stat else 'Outcentro'])
                            ful_unit_bed.write(line+'\n')

                        else:
                            unit_name='nful_'+ltr_name
                            unit_length=sum([int(x[list(x)[0]][2])-int(x[list(x)[0]][1]) for x in target_motif[c_index:c_index+len(unit_list)]])
                            motif_chr_dic=append_dic(chr,unit_name,unit_length,motif_chr_dic)
                        index+=len(unit_list)
                    else:
                        unit_name='nful_'+ltr_name
                        unit_length=sum([int(x[list(x)[0]][2])-int(x[list(x)[0]][1]) for x in target_motif[c_index:c_index+1]])
                        motif_chr_dic=append_dic(chr,unit_name,unit_length,motif_chr_dic)
                        index+=1
                else:
                    unit_name='nful_'+ltr_name
                    unit_length=sum([int(x[list(x)[0]][2])-int(x[list(x)[0]][1]) for x in target_motif[c_index:c_index+1]])
                    motif_chr_dic=append_dic(chr,unit_name,unit_length,motif_chr_dic)
                    # print(unit_name,unit_gap)
                    index+=1
            else:
                # motif_chr_dic=append_dic(chr,seq,length,motif_chr_dic)
                index+=1

        print('-'*20,'stat concat','-'*20)
        Top_5_sp=pd.DataFrame()
        for i in motif_chr_dic:
            chr_dic={}
            top_keys = []
            top_freq = []
            count_dict = {key: value for key, value in motif_chr_dic[i].items() if key.startswith('ful')}
            sorted_keys = sorted(count_dict, key=lambda x: count_dict[x][0], reverse=True)[:5]
            top_freq.extend([count_dict[x] for x in sorted_keys])
            top_keys.extend(sorted_keys)
            num=0
            for x in top_keys:
                num+=1
                chr_dic[i+'_'+str(num)]=x+'__'+str(motif_chr_dic[i][x][0])+'/'+str(motif_chr_dic[i][x][1])
            keys_line=pd.Series(chr_dic,dtype='str')
            Top_5_sp=pd.concat([Top_5_sp,keys_line])
            if fol=='freq':
                line=pd.Series(motif_chr_dic[i]).apply(lambda x: x[0]).to_frame().T
            else:
                line=pd.Series(motif_chr_dic[i]).apply(lambda x: x[1]).to_frame().T
            line.index = pd.Series([i+'_'+sp])
            line['unit_number'] = unit_numb[i]
            final_df=pd.concat([final_df,line])
        for i in intact_motif_chr_dic:
            if fol=='freq':
                line=pd.Series(intact_motif_chr_dic[i]).apply(lambda x: x[0]).to_frame().T
            else:
                line=pd.Series(intact_motif_chr_dic[i]).apply(lambda x: x[1]).to_frame().T
            line.index = pd.Series([i+'_'+sp])
            line['unit_number'] = insert_numb[i]
            intact_df=pd.concat([intact_df,line])
        Top_5_sp=Top_5_sp.rename(columns={0:sp})
        # print(final_df)
        ful_top5_df=pd.concat([ful_top5_df,Top_5_sp],axis=1)
    ful_unit_bed.close()
    ful_top5_df=ful_top5_df.sort_index()
    print(ful_top5_df)
    print('-'*20,'divide line','-'*20)
    final_df.fillna('0',inplace=True)
    # print(final_df.loc[:,'unit_number'])
    final_df=final_df.loc[:,final_df.apply(criterion_freq)]
    intact_df.fillna('0',inplace=True)
    intact_df=intact_df.loc[:,intact_df.apply(criterion_freq2)]
    print(sp,'Chr repeat stat done!')
    print(sp,'Chr intact repeat stat done!')

    if centro_stat:
        if fol=='freq':
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_centro_LTR_freq.csv')
            intact_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_centro_fullLTR_freq.csv')
        else:
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_centro_LTR_length.csv')
            intact_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_centro_fullLTR_length.csv')
    else:
        if fol=='freq':
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_OUTcentro_LTR_freq.csv')
            intact_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_OUTcentro_fullLTR_freq.csv')
        else:
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_OUTcentro_LTR_length.csv')
            intact_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\0123_statresult\01_OUTcentro_fullLTR_length.csv')

    # outfile.close()
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff_rmfup.list"
    # pan centro bed
    # Chr01_13-65 16749781 17724890
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed"
    # pan_centro_bed_file=sys.argv[2]
    # inter_v_ltr=sys.argv[4]
    inter_v_ltr=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\intact_seq_id.txt"
    # ltr_length_file=sys.argv[5]
    ltr_length_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\RepeatMask_16_pan.lib.len"
    # centro_state=sys.argv[3]
    centro_state=0
    if int(centro_state) == 0:
        centro_stat=False
    else:
        centro_stat=True
    # fol in 'freq' or 'length'
    # fol=sys.argv[6]
    fol='length'
    main(pan_gff_file_list,pan_centro_bed_file,True,inter_v_ltr,ltr_length_file,fol)
    main(pan_gff_file_list,pan_centro_bed_file,False,inter_v_ltr,ltr_length_file,fol)
    fol='freq'
    main(pan_gff_file_list,pan_centro_bed_file,True,inter_v_ltr,ltr_length_file,fol)
    main(pan_gff_file_list,pan_centro_bed_file,False,inter_v_ltr,ltr_length_file,fol)