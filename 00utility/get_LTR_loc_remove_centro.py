'''
Descripttion: 修改gff3文件中的loc信息,以centro区域起点为1,即将注释record的起点从序列，更改为CR区域起点
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-03-17 19:39:14
'''
import sys
import re
import logging

def main(LTRfile,centrofile,outfile='LTR_minusCR.bed'):
    pattern=re.compile(r'\"Motif:(.*)\"')
    outfile=open(outfile,'w')
    centro_loc={}
    with open(centrofile,'r') as mf:
        for i in mf:
            if 'start' in i:
                continue
            line=i.strip().split()
            centro_loc[line[0]]=[int(line[1]),int(line[2])]
    with open(LTRfile,'r') as LTRf:
        for i in LTRf:
            line=i.strip().split()
            chr=line[0]
            centro_s=centro_loc[chr][0]
            centro_e=centro_loc[chr][1]
            seq_name=pattern.findall(line[7])[0]
            if int(line[1]) > centro_e or int(line[2]) < centro_s or 'SAT-2_OS' in i:
                continue
            else:
                line[0]=line[0]+'_centro'
                line[1]=str(int(line[1])-centro_s)
                line[2]=str(int(line[2])-centro_s)
                line[7]=seq_name
                write_line='\t'.join(line[:3]+line[7:8]+line[3:5]+line[8:10])
            # print(write_line)
            outfile.write(write_line+'\n')
    outfile.close()

if __name__ == "__main__":
    LTR_gff_file=sys.argv[1]
    #LTR_gff_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\merge_full_LTR.bed'
    centro_loc_file=sys.argv[2]
    #centro_loc_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    outfile=sys.argv[3]
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(LTR_gff_file,centro_loc_file,outfile)
