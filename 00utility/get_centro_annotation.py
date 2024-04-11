'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-23 11:08:21
LastEditors: Ne0tea
LastEditTime: 2023-12-19 13:22:37
'''
import re

def main(file,pan_centro_bed_file):
    outfile=open(file.split('.')[0]+'_centro.gff3','w')
    chr_centro_bed_dic={}
    with open(pan_centro_bed_file,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            if line[0].split('_')[0] in chr_centro_bed_dic:
                chr_centro_bed_dic[line[0].split('_')[0]][line[0].split('_')[1]]=[int(line[1]),int(line[2])]
            else:
                chr_centro_bed_dic[line[0].split('_')[0]]={line[0].split('_')[1]:[int(line[1]),int(line[2])]}
    pattern=re.compile(r'\"Motif:(.*)\"')
    ful_seq_name=[]
    with open(file,'r') as gffile:
        for i in gffile:
            # print(i)
            if i.startswith('#'):
                continue
            line=i.strip().split()
            chr=line[0].split('_')[0]
            sp=line[0].split('_')[1]
            if (chr=='Chr07' and sp=='CW06') or int(line[3]) > chr_centro_bed_dic[chr][sp][1] \
                or int(line[4]) < chr_centro_bed_dic[chr][sp][0] :
                continue
            try:
                seq_name=pattern.findall(line[9])[0]
                ful_seq_name.append(seq_name)
            except IndexError:
                # print(line)
                print(line[9])
            outfile.write(i)
    outfile.close()
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pan_gff_file_list=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\pan_gff.list"
    # pan centro bed
    # Chr01_13-65 16749781 17724890
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed"
    with open(pan_gff_file_list,'r') as gfflist:
        for i in gfflist:
            file=i.strip()
            main(file,pan_centro_bed_file)