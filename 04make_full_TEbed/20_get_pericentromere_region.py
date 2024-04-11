'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-30 14:50:52
LastEditors: Ne0tea
LastEditTime: 2024-02-02 23:20:13
'''
def main(orgin_centro_file,fasta_fai,outfile):
    chr_length_dic={}
    with open(fasta_fai,'r') as faif:
        for i in faif:
            line=i.strip().split()
            chr_length_dic[line[0]]=int(line[1])
    # print(chr_length_dic)
    of=open(outfile,'w')
    with open(orgin_centro_file,'r') as cf:
        for i in cf:
            line=i.strip().split()
            cen_start= int(line[1])
            cen_end=int(line[2])
            left_region=chr_length_dic[line[0]]-cen_end
            new_start=int(cen_start-10*cen_start/120)
            new_end=int(cen_end+10*left_region/120)
            out_line=line[0]+'\t'+str(new_start)+'\t'+str(new_end)
            # print(out_line)
            of.write(out_line+'\n')
    of.close()

if __name__ == "__main__":
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed"
    fasta_fai=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\all_asm_70.fa.fai'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_periregion_70m.bed'
    main(pan_centro_bed_file,fasta_fai,outfile)