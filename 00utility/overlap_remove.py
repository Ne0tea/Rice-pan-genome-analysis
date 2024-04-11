'''
Descripttion: V1:针对RM结果中的GFF3文件。合并注释结果中，因为gap被注释成两个相同LTR；去除在相同位置注释成不同LTR的record(取其中置信度最高，长度最长的注释)。
              输出去除这些record的gff3.(rmfup_replaced.gff3)
              V2:合并重复条目中的LTR区域
Author: Ne0tea
version: V2
Date: 2023-11-27 00:09:51
LastEditors: Ne0tea
LastEditTime: 2024-01-24 13:21:54
'''
import re
import os
import linecache
from  multiprocessing import Process

def produce_dedup_gff3(file,seq_len,pattern,sp):
    with open(file,'r') as gffile:
        last_start=0
        last_end=0
        last_seqname=''
        last_qual=100
        last_ali_len=0
        last_ali_start=0
        last_ali_end=0
        write_line=[]
        overlap_list=[]
        split_list=[]
        count=0
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
            if seq_name not in seq_len:
                seq_concus_len=end-start
            else:
                seq_concus_len=seq_len[seq_name]
            if seq_name == last_seqname:
                # print(seq_concus_len)
                if (ali_start > last_ali_end - 0.1*seq_concus_len) or \
                    (ali_end < last_ali_start + 0.1*seq_concus_len):
                    if count-1 in write_line:
                        write_line.remove(count-1)
                        split_list.append(count-1)
                write_line.append(count)
            else:
                if start < last_end - 0.1*seq_concus_len and count>0:
                    if qual < last_qual and ali_len > last_ali_len:
                        # print(count)
                        if count-1 in write_line:
                            write_line.remove(count-1)
                            overlap_list.append(count-1)
                        else:
                            pass
                        write_line.append(count)
                    else:
                        overlap_list.append(count)
                else:
                    write_line.append(count)

            last_start=start
            last_end=end
            last_qual=qual
            last_ali_start=ali_start
            last_ali_end=ali_end
            last_ali_len=ali_len
            last_seqname=seq_name
            count+=1
    outfile=open(file.split('.')[0]+'rmfup_replaced.gff3','w')
    print(sp,'stat')
    for ll in write_line:
        merge_count=[x for x in split_list if x < ll]
        start_list=[]
        end_list=[]
        for i in merge_count:
            i_line=linecache.getline(file, i+1)
            start_list.append(int(i_line.strip().split()[3]))
            end_list.append(int(i_line.strip().split()[4]))
            split_list.remove(i)
        ll_line = linecache.getline(file, ll+1)
        start_list.append(int(ll_line.strip().split()[3]))
        end_list.append(int(ll_line.strip().split()[4]))

        cur_start,cur_end=min(start_list),max(end_list)
        ll_line=ll_line.strip().split()
        ll_line[3:5]=[str(cur_start),str(cur_end)]
        outfile.write('\t'.join(ll_line)+'\n')
    print(sp,'done')
    print('--------------------------------------------------------')
    outfile.close()


if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff.list"
    ltr_length_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\RepeatMask_16_pan.lib.len"
    gff_list=[]
    with open(pan_gff_file_list,'r') as gf:
        for i in gf:
            line=i.strip()
            gff_list.append(line)
    seq_len={}
    with open(ltr_length_file,'r') as lenf:
        for i in lenf:
            line=i.strip().split()
            seq_len[line[0]]=int(line[1])
    pattern=re.compile(r'\"Motif:(.*)\"')
    p_list=[]
    for file in gff_list:
        sp=os.path.basename(file).split('.')[0]
        print(sp,'Start')
        p = Process(target=produce_dedup_gff3, args=(file,seq_len,pattern,sp,))
        p_list.append(p)
        p.start()
    for i in p_list:
        p.join()