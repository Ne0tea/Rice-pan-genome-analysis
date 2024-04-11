'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-23 17:05:40
LastEditors: Ne0tea
LastEditTime: 2023-11-30 19:57:42
'''
import pandas as pd
pd.set_option('display.max_rows', None)

def main(seq_conut_csv,pan_centro_bed_file):
    chr_centro_length={}
    chr_centro_bed_dic={}
    with open(pan_centro_bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            if line[0].split('_')[0] in chr_centro_bed_dic:
                chr_centro_bed_dic[line[0].split('_')[0]][line[0].split('_')[1]]=[int(line[1]),int(line[2])]
                chr_centro_length[line[0].split('_')[0]][line[0].split('_')[1]]=int(line[2])-int(line[1])
            else:
                chr_centro_bed_dic[line[0].split('_')[0]]={line[0].split('_')[1]:[int(line[1]),int(line[2])]}
                chr_centro_length[line[0].split('_')[0]]={line[0].split('_')[1]:int(line[2])-int(line[1])}
    seq_stat_df=pd.read_csv(seq_conut_csv,index_col=0,header=0)
    # print(seq_stat_df)
    for rownames,item in seq_stat_df.iterrows():
        frag_item=item[[x for x in item.index.to_list() if x.split('_')[0] =='frag']].str.split('/',expand=True)[1].to_frame().astype('int')
        ful_item=item[[x for x in item.index.to_list() if x.split('_')[0] =='ful']].str.split('/',expand=True)[1].astype('int')
        insert_item=item[[x for x in item.index.to_list() if x.split('_')[0] =='insert']].str.split('/',expand=True)[1].to_frame().astype('int')
        # line=item.str.split('/',expand=True)[1].to_frame().astype('int')
        print(rownames,chr_centro_length[rownames.split('_')[0]][rownames.split('_')[1]],
                ful_item.nlargest(5).index.tolist(),
                ful_item.nlargest(5).tolist())
        # print(rownames,insert_item.idxmax(axis=0))
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    seq_conut_csv=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\1130_statresult\Chr_seq_centro_region_intact.csv"
    # pan centro bed
    # Chr01_13-65 16749781 17724890
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\70_pan_centro.bed"
    # pan_centro_bed_file=sys.argv[2]
    main(seq_conut_csv,pan_centro_bed_file)