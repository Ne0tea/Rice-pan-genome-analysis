'''
Descripttion: get target superfamily annotation records from gff3
Author: Ne0tea
version: 
Date: 2023-11-26 16:26:50
LastEditors: Ne0tea
LastEditTime: 2023-12-21 11:32:10
'''
import pandas as pd
import re
import pandas as pd
import os
pd.set_option('display.max_columns', None)

def criterion_length(series):
    state=False
    if series.name == 'unit_number':
        return True
    for i in series:
        if int(i) >= 10000 :
            state=True
            return state
    return state
def main(gff3_file,pan_centro_bed_file,repeat_id,target_family,sp_phy_file):
    gfflist=[]
    with open(gff3_file,'r') as gf:
        for i in gf:
            gfflist.append(i.strip())
    pattern=re.compile(r'\"Motif:(.*)\"')
    # outfile=open(gff3_file.split('.')[0]+'_target'+target_family+'.bed','w')
    sp_phy={}
    with open(sp_phy_file,'r') as phyf:
        for i in phyf:
            line=i.strip().split()
            sp_phy[line[0]]=line[1]
    sp_centro_length_dic={}
    with open(pan_centro_bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            if line[0].split('_')[1] in sp_centro_length_dic:
                sp_centro_length_dic[line[0].split('_')[1]]+=int(line[2])-int(line[1])
            else:
                sp_centro_length_dic[line[0].split('_')[1]]=int(line[2])-int(line[1])
    type_dic={}
    with open(repeat_id,'r') as idf:
        for i in idf:
            # print(pattern.findall(i.strip().split('#')[0]))
            if '#' not in i:
                continue
            # print(i.strip())
            seq_name=re.sub(r'_OS|-I|-LTR|_LTR|_I','',i.strip().split('#')[0])
            # print(i.strip().split('#'))
            class1=i.strip().split('#')[1].split('/')[0]
            try:
                class2=i.strip().split('#')[1].split('/')[1]
            except IndexError:
                # print(seq_name)
                class2='None'
            type_dic[seq_name]=[class1,class2]
    # print(type_dic)
    final_sp_df=pd.DataFrame()
    for gfffile in gfflist:
        seq_length_dic={}
        seq_order_dic={}
        with open(gfffile,'r') as af:
            for i in af:
                if i.startswith('#'):
                    continue
                line=i.strip().split()
                sp=line[0].split('_')[1]
                try:
                    seq_name=re.sub(r'_OS|-I|-LTR|_LTR|_I','',pattern.findall(line[9])[0])
                except IndexError:
                    print(seq_name)
                    # continue
                    pass
                if seq_name in type_dic and type_dic[seq_name][1]==target_family:
                    length=int(line[4])-int(line[3])
                    if seq_name in seq_order_dic:
                        seq_order_dic[seq_name]+=1
                    else:
                        seq_order_dic[seq_name]=0
                    # outline=line[0]+'\t'+line[3]+'\t'+line[4]+'\t'+sp+'_'+seq_name+'_'+str(seq_order_dic[seq_name])+\
                    #         '\t'+'0'+'\t'+line[6]+'\t'+pattern.findall(line[9])[0]+'\t'+sp+'\t'+sp_phy[sp]+'\n'
                    # outfile.write(outline)
                    if seq_name in seq_length_dic:
                        seq_length_dic[seq_name]+=length
                    else:
                        seq_length_dic[seq_name]=length
        cur_df=pd.Series(seq_length_dic).to_frame().T
        # print(sp)
        cur_df.index=pd.Series([sp])
        final_sp_df=pd.concat([final_sp_df,cur_df])
    final_sp_df.fillna(0,inplace=True)
    # final_sp_df=final_sp_df.loc[:,final_sp_df.apply(criterion_length)]
    final_sp_df = final_sp_df.astype(int)
    Other_line=pd.Series(sp_centro_length_dic).to_frame()
    Other_line.columns=pd.Series(['centro_length'])
    # print(final_sp_df.filter(regex='SZ').loc[:, ~final_sp_df.filter(regex='SZ').columns.str.contains('SZ-22')])
    
    final_sp_df['RIRE_family'] = final_sp_df.filter(regex='RIRE').sum(axis=1,numeric_only=True)
    final_sp_df['Gypsy_family'] = final_sp_df.filter(regex='Gypsy').sum(axis=1,numeric_only=True)
    final_sp_df['ATL_family'] = final_sp_df.filter(regex='ATL').sum(axis=1,numeric_only=True)
    final_sp_df['SZ_family'] = final_sp_df.filter(regex='SZ').loc[:, ~final_sp_df.filter(regex='SZ').columns.str.contains('SZ-22')].sum(axis=1,numeric_only=True)
    final_sp_df['RETRO_family'] = final_sp_df.filter(regex='RETRO').sum(axis=1,numeric_only=True)
    final_sp_df['TRUNCATOR_family'] = final_sp_df.filter(regex='TRUNCATOR').sum(axis=1,numeric_only=True)
    final_sp_df=final_sp_df.filter(regex='family|SZ-22').loc[:,['SZ-22','RIRE_family','Gypsy_family','ATL_family','TRUNCATOR_family','SZ_family','RETRO_family']]
    # final_sp_df['Others']=Other_line.loc[final_sp_df.index].sum(axis=1)-final_sp_df.sum(axis=1,numeric_only=True)
    final_sp_df.insert(0,'loc',-1)
    final_sp_df.insert(1,'radius',final_sp_df.iloc[:, 1:].sum(axis=1)/200000)
    final_sp_df.to_csv(r'E:\Bio_analysis\Weedyrice\pan_weedyrice\1207_statresult\Sp_Seqname_ratio.txt',sep='\t')

    # print(final_sp_df)
    # outfile.close()
if __name__ == "__main__":

    # gff3_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\Chr11_centro.gff3"
    centro_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_centro_gff.list"
    sp_phy_file=r'pan_weedyrice/70_phylogenic.txt'
    repeat_id=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\repeatmask_lib.id'
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed"
    target_family='Gypsy'
    main(centro_file_list,pan_centro_bed_file,repeat_id,target_family,sp_phy_file)
    # main(gff3_file,repeat_id,target_family,sp_phy_file)