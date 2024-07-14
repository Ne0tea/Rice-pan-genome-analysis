'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-02-13 12:44:22
'''
import sys
import logging
import re
import os
import pandas as pd
# pd.set_option('display.max_rows',None)
def contain_info(row):
    strings=row['name']
    start_site=int(row['queryStart'])
    end_site=int(row['queryEnd'])
    # print(strings,start_site,end_site)
    if strings.endswith('fw'):
        if end_site > 225:
            return True
    elif strings.endswith('aw'):
        if start_site < 75:
            return True
    return False

def calculate_insert_loc(df):
    ref_start=df['refStart'].to_list()
    ref_end=df['refEnd'].to_list()
    if (ref_start[0] < ref_end[0] and ref_start[1] < ref_end[1]):
        insert_loc=(max(ref_start)+min(ref_end))/2
        return 'Common',insert_loc
    elif ref_start[0] > ref_end[0] and ref_start[1] > ref_end[1]:
        insert_loc=(max(ref_end)+min(ref_start))/2
        return 'Common',insert_loc
    else:
        loc_line = df['identity'].max()
        loc_row=df[df['identity'] == loc_line]
        #print(loc_row.iloc[0][7])
        if loc_row['name'].values[0].endswith('fw'):
            insert_loc=int(loc_row.iloc[0][7])
        elif loc_row['name'].values[0].endswith('aw'):
            insert_loc=int(loc_row.iloc[0][6])
        return 'InvertRepeat',insert_loc

def main(blast_file):
    Invert_number=0
    total_number=0
    fw_dataframe=pd.DataFrame(columns=['name','monomer','identity','aligmentlength','queryStart','queryEnd','refStart','refEnd'])
    aw_dataframe=pd.DataFrame(columns=['name','monomer','identity','aligmentlength','queryStart','queryEnd','refStart','refEnd'])
    result_dataframe=pd.DataFrame(columns=['name','monomer','identity','aligmentlength','queryStart','queryEnd','refStart','refEnd'])
    insert_loc_outf=open('SZ-22_insert_loc.list','w')
    ID_monomer_outf=open('SZ-22_ID_monomer.txt','w')
    ful_dataframe=pd.read_table(blast_file,usecols=[0,1,2,3,6,7,8,9],names=['name','monomer','identity','aligmentlength','queryStart','queryEnd','refStart','refEnd'],sep='\t')
    filter_dataframe=ful_dataframe[ful_dataframe[['name','queryStart','queryEnd']].apply(lambda x : contain_info(x),axis=1)]
    fw_dataframe=filter_dataframe.query('name.str.endswith("fw")')
    aw_dataframe=filter_dataframe.query('name.str.endswith("aw")')
    for i in fw_dataframe['name'].unique():
        aw_name=i.replace('fw','aw')
        fw_cur_monomer_id=fw_dataframe.loc[fw_dataframe['name']==i,'monomer'].to_list()
        aw_cur_monomer_id=aw_dataframe.loc[aw_dataframe['name']==aw_name,'monomer'].to_list()
        info_monomer=list(set(fw_cur_monomer_id) & set(aw_cur_monomer_id))
        # if i=='Chr01_13-65_SZ-22_LTR5_fw':
        #     print(fw_dataframe.loc[fw_dataframe['name']==i,:])
        #     print(aw_dataframe.loc[fw_dataframe['name']==i,:])
        monomer_identity_dic={}
        fw_idic={}
        aw_idic={}
        unit_dataframe=pd.DataFrame(columns=['name','monomer','identity','aligmentlength','queryStart','queryEnd','refStart','refEnd'])
        # print(info_monomer)
        if info_monomer:
            fw_curdf=fw_dataframe.loc[fw_dataframe['name']==i,:]
            fw_curdf = fw_curdf[fw_curdf['monomer'].isin(info_monomer).to_list()]
            for monomer,ident in zip(fw_curdf['monomer'],fw_curdf['identity']):
                ident=float(ident)
                if monomer in fw_idic:
                    pass
                    fw_idic[monomer]=max(fw_idic[monomer],ident)
                else:
                    fw_idic[monomer]=ident

            aw_curdf=aw_dataframe.loc[aw_dataframe['name']==aw_name,:]
            aw_curdf = aw_curdf[aw_curdf['monomer'].isin(info_monomer).to_list()]

            for monomer,ident in zip(aw_curdf['monomer'],aw_curdf['identity']):
                ident=float(ident)
                if monomer in aw_idic:
                    pass
                    aw_idic[monomer]=max(aw_idic[monomer],ident)
                else:
                    aw_idic[monomer]=ident

            for x in fw_idic:
                monomer_identity_dic[x]=fw_idic[x]+aw_idic[x]
            cur_monomer_name=max(monomer_identity_dic,key=lambda x:monomer_identity_dic[x])
            fw_max_row=fw_curdf.loc[fw_curdf['monomer']==cur_monomer_name,:]

            recent_date = fw_max_row['queryEnd'].max()
            fw_max_row=fw_max_row[fw_max_row['queryEnd'] == recent_date]

            aw_max_row=aw_curdf.loc[aw_curdf['monomer']==cur_monomer_name,:]
            recent_date = fw_max_row['queryStart'].min()
            fw_max_row=fw_max_row[fw_max_row['queryStart'] == recent_date]

            unit_dataframe=pd.concat([fw_max_row,aw_max_row],axis=0)
            result_dataframe=pd.concat([result_dataframe,unit_dataframe],axis=0)
            status,insert_loc=calculate_insert_loc(unit_dataframe)
            if status=='InvertRepeat':
                Invert_number+=1
            total_number+=1
            insert_loc_outf.write(str(insert_loc)+'\n')
            print(i,'calculate complete')
            monomer_line=i+'\t'+cur_monomer_name+'\t'+str(insert_loc)
            ID_monomer_outf.write(monomer_line+'\n')
    insert_loc_outf.close()
    ID_monomer_outf.close()
    print(Invert_number,total_number)
if __name__ == "__main__":
    # monomer_file=sys.argv[1]#SZ-22_unit_centro.bed 
    # monomer_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\monomer_anno.txt'
    # blast_file=sys.argv[1]
    blast_file=r'D:\study resource\13-65_SZ-22.blast_out'

    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(blast_file)
