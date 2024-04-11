'''
Descripttion:针对RM结果中的GFF3文件。合并注释结果中，因为gap被注释成两个相同LTR(允许0.1concensus length的幅度)；去除在相同位置注释成不同LTR的record(取其中置信度最高，长度最长的注释)。
             输出在这些record背后添加标签的gff3.(rmfup.gff3)
Author: Ne0tea
version: 
Date: 2023-11-27 00:09:51
LastEditors: Ne0tea
LastEditTime: 2023-12-01 14:22:47
'''
import re
import os
from collections import Counter
import logging

def main(pan_gff_file_list,ltr_len_file):
    gff_list=[]
    with open(pan_gff_file_list,'r') as gf:
        for i in gf:
            line=i.strip()
            gff_list.append(line)
    seq_len={}
    with open(ltr_len_file,'r') as lenf:
        for i in lenf:
            line=i.strip().split()
            seq_len[line[0]]=int(line[1])
    pattern=re.compile(r'\"Motif:(.*)\"')
    for file in gff_list:
        sp=os.path.basename(file).split('.')[0]
        print(sp,'Start')
        logging.info('%s Stat Start',sp)
        logging.info('-'*20)
        with open(file,'r') as gffile:
            last_start=0
            last_end=0
            last_seqname=''
            last_qual=100
            last_ali_len=0
            last_ali_start=0
            last_ali_end=0
            write_line=[]
            count=0
            homo_anno_seq=[]
            outfile=open(file.split('.')[0]+'rmfup.gff3','w')
            for i in gffile:
                if i.startswith('#'):
                    count+=1
                    continue
                line=i.strip().split()
                '''
                length=int(line[4])-int(line[3])
                '''
                start=int(line[3])
                end=int(line[4])
                chr_name=line[0].split('_')[0]
                ali_start=int(line[10])
                ali_end=int(line[11])
                ali_len=int(line[11])-int(line[10])
                qual=float(line[5])
                seq_name=pattern.findall(line[9])[0]
                if count%100000==0:
                    print(str(count/10000)+'W record passed')
                if seq_name not in seq_len:
                    seq_concus_len=end-start
                else:
                    seq_concus_len=seq_len[seq_name]
                if seq_name == last_seqname:
                    if (ali_start > last_ali_end - 0.1*seq_concus_len) or \
                        (ali_end < last_ali_start + 0.1*seq_concus_len):
                        if count-1 in write_line:
                            write_line.remove(count-1)
                            outfile.write('\tMerged')
                    write_line.append(count)
                    outfile.write('\n')
                    outfile.write(i.strip())
                else:
                    if start < last_end - 0.1*seq_concus_len and count>0:
                        if qual < last_qual and ali_len > last_ali_len:
                            # print(count)
                            if count-1 in write_line:
                                write_line.remove(count-1)
                                homo_anno_seq.append(last_seqname)
                                outfile.write('\tMisanno\n')
                            else:
                                outfile.write('\n')
                                pass
                            write_line.append(count)
                            outfile.write(i.strip())
                        else:
                            if count>2:
                                outfile.write('\n')
                            outfile.write(i.strip()+'\tMisanno')
                            homo_anno_seq.append(seq_name)
                    else:
                        write_line.append(count)
                        if count>2:
                            outfile.write('\n')
                        outfile.write(i.strip())
                last_start=start
                last_end=end
                last_qual=qual
                last_ali_start=ali_start
                last_ali_end=ali_end
                last_ali_len=ali_len
                last_seqname=seq_name
                count+=1
            outfile.close()
        logging.info(Counter(homo_anno_seq).most_common(20))
        logging.info('-'*20)
        print(sp,'done')
        print('-'*20)
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff.list"
    ltr_length_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\RepeatMask_16_pan.lib.len"

    logging.basicConfig(filename='70_pan_remove_overlap.log', \
                        level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s',\
                        filemode='w')
    main(pan_gff_file_list,ltr_length_file)