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
Date: 2023-11-16 19:46:32
LastEditors: Ne0tea
LastEditTime: 2023-12-18 22:47:59
'''
import sys
import re
import pandas as pd
import os
import logging
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def criterion_freq(series):
    state=False
    global fol
    # print(series)
    if series.name == 'LIL_number':
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
    if series.name == 'LIL_number':
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

def find_intersecting_ranges(ranges):
    events = []
    for i, range_item in enumerate(ranges):
        start, end = range_item
        events.append((start, 1, i))
        events.append((end, -1, i))
    events.sort(key=lambda x: (x[0], -x[1], x[2]))
    count = 0
    intersecting_indices = set()
    for _, event_type, index in events:
        if event_type == 1:
            count += 1
        else:
            count -= 1
        if count > 1:
            intersecting_indices.add(index)
    return list(intersecting_indices)

def main(pan_gff_file_list,pan_centro_bed_file,centro_stat,inter_v_ltr:str,ltr_length_file,fol):
    # outfile=open('centro_region_seq.gff3','w')
    ivl_dic={}
    final_df=pd.DataFrame()
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
    for file in gff_list:
        target_motif_chr={}
        sp=os.path.basename(file).split('.')[0].replace('rmfup_replaced','')
        motif_chr_dic={}#{seq_name:[chr_name,start,end,gap_ratio]}
        ltr_numb={}
        with open(file,'r') as gffile:
            count=0
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
                if centro_stat:
                    if (chr_name=='Chr07' and sp=='CW06') or int(line[3]) > chr_centro_bed_dic[chr_name][sp][1] \
                        or int(line[4]) < chr_centro_bed_dic[chr_name][sp][0] :
                        continue
                else:
                    # print(1)
                    if (chr_name=='Chr07' and sp=='CW06') or (int(line[4]) < chr_centro_bed_dic[chr_name][sp][1] and int(line[4]) > chr_centro_bed_dic[chr_name][sp][0]) or \
                        (int(line[3]) > chr_centro_bed_dic[chr_name][sp][0] and int(line[3]) < chr_centro_bed_dic[chr_name][sp][1]) :
                        continue
                # gap_ratio=round(last_gap/seq_concus_len,3)
                if chr_name in target_motif_chr:
                    target_motif_chr[chr_name].append({seq_name:[start,end]})
                else:
                    target_motif_chr[chr_name]=[{seq_name:[start,end]}]
                count+=1
                last_loc=end
        logging.info('%s running',sp)
        for chr in target_motif_chr:
            if chr not in ltr_numb:
                ltr_numb[chr]=0
            index=0
            span=200
            while index < len(target_motif_chr[chr])-1:
                seq=list(target_motif_chr[chr][index])[0]
                start=int(target_motif_chr[chr][index][seq][0])
                end=int(target_motif_chr[chr][index][seq][1])
                length=int(target_motif_chr[chr][index][seq][1])-int(target_motif_chr[chr][index][seq][0])
                # if index % 10000==0:
                #     logging.info('%d0W solo count done',index/10000)
                # print(seq)
                if seq in ivl_dic:
                    c_index=index
                    # print('in ivl',c_index)
                    if 'LTR' in seq:
                        ltr_name=seq
                        inter_name=ivl_dic[seq]
                    else:
                        ful_state=False
                        inter_name=seq
                        ltr_name=ivl_dic[seq]
                    if c_index > 5 and c_index < len(target_motif_chr)-5:
                        cur_ltr_list=target_motif_chr[chr][c_index-6:c_index+5]#[{'A':[0,199]},{'B':[200,300]},{'C':[500,600]}]
                        cur_ltr_range=[cur_ltr_list[0][list(cur_ltr_list[0])[0]][0],cur_ltr_list[-1][list(cur_ltr_list[-1])[0]][1]]
                        left_l=start-span if start > span else 0
                        right_l=end+span 
                        if cur_ltr_range[1]-cur_ltr_range[0] < span:
                            # print(left_l,right_l,cur_ltr_range)
                            print(left_l,right_l,cur_ltr_range,'Error: there is not enough')
                        cur_ltr_loc_list=[x[list(x)[0]] for x in cur_ltr_list]
                        cur_ltr_loc_list[6]=[left_l,right_l]
                        overlap_seq=[list(cur_ltr_list[i])[0] for i in find_intersecting_ranges(cur_ltr_loc_list)]
                        if overlap_seq.count(ltr_name) + overlap_seq.count(inter_name) > 1:
                            unit_name='nsolo_'+ltr_name
                            motif_chr_dic=append_dic(chr,unit_name,length,motif_chr_dic)
                            index+=1
                        else:
                            unit_name='solo_'+ltr_name
                            motif_chr_dic=append_dic(chr,unit_name,length,motif_chr_dic)
                            index+=1
                    elif c_index <= 5:
                        cur_ltr_list=target_motif_chr[chr][0:c_index+5]#[{'A':[0,199]},{'B':[200,300]},{'C':[500,600]}]
                        cur_ltr_range=[cur_ltr_list[0][list(cur_ltr_list[0])[0]][0],cur_ltr_list[-1][list(cur_ltr_list[-1])[0]][1]]
                        left_l=start-span if c_index != 0 else start
                        right_l=end+span 
                        if cur_ltr_range[1]-cur_ltr_range[0] < span:
                            # print(left_l,right_l,cur_ltr_range)
                            # print(cur_ltr_list)
                            print(left_l,right_l,cur_ltr_range,'Error: there is not enough')
                        cur_ltr_loc_list=[x[list(x)[0]] for x in cur_ltr_list]
                        cur_ltr_loc_list[c_index]=[start-span,end+span]
                        overlap_seq=[list(cur_ltr_list[i])[0] for i in find_intersecting_ranges(cur_ltr_loc_list)]
                        if overlap_seq.count(ltr_name) + overlap_seq.count(inter_name) > 1:
                            unit_name='nsolo_'+ltr_name
                            motif_chr_dic=append_dic(chr,unit_name,length,motif_chr_dic)
                            index+=1
                        else:
                            unit_name='solo_'+ltr_name
                            motif_chr_dic=append_dic(chr,unit_name,length,motif_chr_dic)
                            index+=1
                    elif c_index >= len(target_motif_chr)-5:
                        cur_ltr_list=target_motif_chr[chr][c_index-6:]#[{'A':[0,199]},{'B':[200,300]},{'C':[500,600]}]
                        # print(cur_ltr_list)
                        cur_ltr_range=[cur_ltr_list[0][list(cur_ltr_list[0])[0]][0],cur_ltr_list[-1][list(cur_ltr_list[-1])[0]][1]]
                        left_l=start-span
                        right_l=end+span if c_index != len(target_motif_chr)-1 else end
                        if cur_ltr_range[1]-cur_ltr_range[0] < span:
                            # print(left_l,right_l,cur_ltr_range)
                            # print(cur_ltr_list)
                            print(left_l,right_l,cur_ltr_range,'Error: there is not enough')
                        cur_ltr_loc_list=[x[list(x)[0]] for x in cur_ltr_list]
                        # print(cur_ltr_loc_list)
                        cur_ltr_loc_list[6]=[left_l,right_l]
                        overlap_seq=[list(cur_ltr_list[i])[0] for i in find_intersecting_ranges(cur_ltr_loc_list)]
                        if overlap_seq.count(ltr_name) + overlap_seq.count(inter_name) > 1:
                            unit_name='nsolo_'+ltr_name
                            motif_chr_dic=append_dic(chr,unit_name,length,motif_chr_dic)
                            index+=1
                        else:
                            unit_name='solo_'+ltr_name
                            motif_chr_dic=append_dic(chr,unit_name,length,motif_chr_dic)
                            index+=1
                    ltr_numb[chr]+=1
                else:
                    # motif_chr_dic=append_dic(chr,seq,length,motif_chr_dic)
                    index+=1
            # logging.info('%s solo count done',sp)
        print('-'*20,'stat concat','-'*20)
        Top_5_sp=pd.DataFrame()
        for i in motif_chr_dic:
            chr_dic={}
            top_keys = []
            top_freq = []
            count_dict = {key: value for key, value in motif_chr_dic[i].items() if key.startswith('solo')}
            sorted_keys = sorted(count_dict, key=lambda x: count_dict[x][0], reverse=True)[:5]
            top_freq.extend([count_dict[x] for x in sorted_keys])
            top_keys.extend(sorted_keys)
            num=0
            for x in top_keys:
                num+=1
                chr_dic[i+'_'+str(num)]=x+'__'+str(motif_chr_dic[i][x][0])+'/'+str(motif_chr_dic[i][x][1])
            keys_line=pd.Series(chr_dic,dtype='str')
            Top_5_sp=pd.concat([Top_5_sp,keys_line])
            # print(motif_chr_dic[i])
            if fol=='freq':
                line=pd.Series(motif_chr_dic[i]).apply(lambda x: x[0]).to_frame().T
            else:
                line=pd.Series(motif_chr_dic[i]).apply(lambda x: x[1]).to_frame().T
            # print(line)
            line['LIL_number']=ltr_numb[chr]
            # print(line)
            line.index = pd.Series([i+'_'+sp])
            final_df=pd.concat([final_df,line])

        Top_5_sp=Top_5_sp.rename(columns={0:sp})
        # print(Top_5_sp.columns)
        ful_top5_df=pd.concat([ful_top5_df,Top_5_sp],axis=1)
    ful_top5_df=ful_top5_df.sort_index()
    print(ful_top5_df)
    print('-'*20,'divide line','-'*20)
    final_df.fillna('0',inplace=True)
    final_df=final_df.loc[:,final_df.apply(criterion_freq)]
    print(sp,'Chr repeat stat done!')

    if centro_stat:
        if fol=='freq':
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_centro_soloLTR_freq.csv')
        else:
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_centro_soloLTR_length.csv')
    else:
        if fol=='freq':
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_OUTcentro_soloLTR_freq.csv')
        else:
            final_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\01_OUTcentro_soloLTR_length.csv')

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
    centro_state=1
    if int(centro_state) == 0:
        centro_stat=False
    else:
        centro_stat=True
    # fol in 'freq' or 'length'
    # fol=sys.argv[6]
    fol='freq'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(pan_gff_file_list,pan_centro_bed_file,True,inter_v_ltr,ltr_length_file,fol)
    main(pan_gff_file_list,pan_centro_bed_file,False,inter_v_ltr,ltr_length_file,fol)