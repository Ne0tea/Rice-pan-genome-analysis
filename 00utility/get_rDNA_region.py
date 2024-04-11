'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-02-26 16:48:03
LastEditors: Ne0tea
LastEditTime: 2024-02-29 21:17:47
'''
import sys
def exists_within_distance(ranges_list, given_range, distance_threshold=500):
    adj_list=[]
    for start, end in ranges_list:
        min_distance = min(abs(given_range[0] - end), abs(start - given_range[1]))
        # 如果最小距离小于等于指定阈值，则存在距离在500bp以内的范围
        if min_distance <= distance_threshold:
            adj_list.append((start,end))
    # 如果列表中不存在符合条件的范围，则返回 False
    if adj_list:
        return adj_list
    else:
        return False

def merge_ranges(ranges):
    # 首先按照范围的起始值进行排序
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    merged_ranges = []
    # 初始化当前范围为第一个范围
    current_range = sorted_ranges[0]
    # 遍历排序后的范围列表
    for i in range(1, len(sorted_ranges)):
        # 当前范围的结束值
        current_end = current_range[1]
        # 下一个范围的起始值和结束值
        next_start, next_end = sorted_ranges[i]
        # 如果当前范围与下一个范围相邻或重叠，则合并两个范围
        if current_end >= next_start - 1:
            current_range = (current_range[0], max(current_end, next_end))
        else:
            # 如果不相邻，则将当前范围加入结果列表，并更新当前范围为下一个范围
            merged_ranges.append(current_range)
            current_range = (next_start, next_end)
    # 将最后一个范围加入结果列表
    merged_ranges.append(current_range)
    return merged_ranges

def main(rDNA_blast_file,rDNA_length):
    outfile=rDNA_blast_file.strip('_blast.out')+'rDNA_region.txt'
    of=open(outfile,'w')
    rDNA_len_dic={}
    with open(rDNA_length,'r') as rf:
        for i in rf:
            line=i.strip().split()
            rDNA_len_dic[line[0]]=int(line[1])
    c_rDNA_dic={}
    with open(rDNA_blast_file,'r') as rbf:
        for i in rbf:
            line=i.strip().split()
            rDNA_name=line[0].split('_')[0]
            rDNA_len=rDNA_len_dic[line[0]]
            if rDNA_len - int(line[3]) > 0.1*rDNA_len or float(line[2]) < 95:
                continue
            ref_start=min(int(line[8]),int(line[9]))
            ref_end=max(int(line[8]),int(line[9]))
            chr_name=line[1]
            if rDNA_name+'_'+chr_name in c_rDNA_dic:
                c_rDNA_dic[rDNA_name+'_'+chr_name].append((ref_start,ref_end))
            else:
                c_rDNA_dic[rDNA_name+'_'+chr_name]=[(ref_start,ref_end)]
    for i in c_rDNA_dic:
        c_rDNA_dic[i]=merge_ranges(c_rDNA_dic[i])
    id_58=[x for x in c_rDNA_dic.keys() if '5.8S' in x]
    id_5=[x for x in c_rDNA_dic.keys() if x.startswith('5S')]
    result_dic={}
    #c_18_id is 18S rDNA locate at one chromosome
    # print(id_58)
    for c_58_id in id_58:
        c_18_id=c_58_id.replace('5.8S','18S')
        c_25_id=c_58_id.replace('5.8S','25S')
        if c_18_id in c_rDNA_dic and c_25_id in c_rDNA_dic:
            c_18=c_rDNA_dic[c_18_id]
            c_58=c_rDNA_dic[c_58_id]
            c_25=c_rDNA_dic[c_25_id]
            chr_name=c_18_id.replace('18S_','')+'_45S_'
            # print(c_18,c_58,c_25)
            count=0
            for i in c_58:
                count+=1
                r45S=[]
                r45S.append(i)
                exist_18=exists_within_distance(c_18, i)
                if not exist_18:
                    continue
                elif len(exist_18) > 1:
                    print(c_18_id.replace('18S_',''),'may have clustered 5.8S result !')
                r45S.extend(exist_18)
                exist_25=exists_within_distance(c_25, i)
                if not exist_25:
                    continue
                elif len(exist_25) > 1:
                    print(c_18_id.replace('18S_',''),'may have clustered 5.8S result !')
                r45S.extend(exist_25)
                r45S=sorted(r45S, key=lambda x: x[0])
                outline=chr_name+str(count)+'\t'+"\t".join([f"{start},{end}" for start, end in r45S])
                # print(outline)
                of.write(outline+'\n')
        else:
            continue

    for i in id_5:
        count=0
        chr_name=i.replace('5S_','')+'_5S_'
        # print(i)
        for x in c_rDNA_dic[i]:
            count+=1
            outline=chr_name+str(count)+'\t'+str(x[0])+','+str(x[1])
            # print(outline)
            of.write(outline+'\n')
if __name__ == "__main__":
    #blast out
    # rDNA_blast_file=sys.argv[1]
    rDNA_blast_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\rDNA\rDNA_MH63_blast.out"
    # rDNA_length=sys.argv[2]
    rDNA_length=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\irgsp_rDNA.fasta.fai'
    main(rDNA_blast_file,rDNA_length)