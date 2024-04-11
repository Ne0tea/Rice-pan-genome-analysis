'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-31 21:30:35
LastEditors: Ne0tea
LastEditTime: 2024-01-31 21:36:58
'''
def main(cluster_anno_file,outfile):
    record_dic={}
    with open(cluster_anno_file,'r') as wf:
        for i in wf:
            line=i.strip().split()
            if line[6] in record_dic:
                old_s=record_dic[line[6]][1]
                old_e=record_dic[line[6]][2]
                new_s=min(record_dic[line[6]][1],line[1])
                new_e=max(record_dic[line[6]][2],line[2])
                record_dic[line[6]][1]=new_s
                record_dic[line[6]][2]=new_e
            else:
                record_dic[line[6]]=[line[0],line[1],line[2],line[5]]
    with open(outfile,'w') as of:
        for i in record_dic:
            line='\t'.join(record_dic[i][:3])+'\t'+i+'\t'+record_dic[i][3]
            of.write(line+'\n')
if __name__ == "__main__":
    cluster_anno_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22_cluster_distribution.bed'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22_cluster_distribution_merge.bed'
    main(cluster_anno_file,outfile)