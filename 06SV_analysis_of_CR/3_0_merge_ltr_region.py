'''
Descripttion: 全基因组ltr注释文件比较大,建议在服务器上运行
USAGE:python3 $0 ltr_anno(should be sorted) centro_file > outfile
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-03-21 20:34:42
'''
import sys
import re
import logging

def cal_range_distance(range1,range2):
    start1,end1=min(range1),max(range1)
    start2,end2=min(range2),max(range2)
    if end1 < start2:
        distance=start2-end1
    elif end2 < start1:
        distance=start1-end2
    else:
        distance=0
    return distance

def main(ltr_file,centro_file):
    centro_dic={}
    with open(centro_file,'r') as cf:
        for i in cf:
            if 'start' in i:
                continue
            line=i.strip().split()
            centro_dic[line[0]]=[int(line[1]),int(line[2])]
    with open(ltr_file,'r') as lf:
        used_chr_sp=[]
        for i in lf:
            line=i.strip().split()
            if 'Motif:SAT' in i:
                continue
            c_start=int(line[1])
            c_end=int(line[2])
            centro_region=centro_dic[line[0]]
            if line[0] not in used_chr_sp:
                if 'culmate_start' in locals() and 'culmate_end' in locals():
                    del culmate_end
                    del culmate_start
                last_end=0
                used_chr_sp.append(line[0])

            if cal_range_distance([c_start,c_end],centro_region) > 2000000:
                last_end=c_end
                continue
            if c_start-last_end > 200:
                if 'culmate_start' in locals() and 'culmate_end' in locals():
                    outline='\t'.join([line[0],str(culmate_start),str(culmate_end),'ltr_region'])
                    print(outline)
                culmate_start=c_start
                culmate_end=c_end
            else:
                culmate_end=c_end
            last_end=c_end
        outline='\t'.join([line[0],str(culmate_start),str(culmate_end),'ltr_region'])
        print(outline)


if __name__ == "__main__":
    ltr_anno_file=sys.argv[1]
    centro_file=sys.argv[2]
    # centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    # ltr_anno_file=r''

    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(ltr_anno_file,centro_file)