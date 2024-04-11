'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-02-26 16:48:03
LastEditors: Ne0tea
LastEditTime: 2024-03-01 19:05:31
'''
def has_intersection(range1, range2):
    range1 = [int(value) for value in range1]
    range2 = [int(value) for value in range2]
    if range1[1] < range2[0] or range2[1] < range1[0]:
        return False
    return True

def main(centro_file,pericentro_LTR_file,outfile):
    centro_region={}
    out=open(outfile,'w')
    with open(centro_file,'r') as cf:
        for i in cf:
            # if i.startswith(''):
            #     continue
            line=i.strip().split()
            centro_region[line[0]]=[int(line[1]),int(line[2])]
    C=0
    OCIP=0
    with open(pericentro_LTR_file,'r') as pf:
        for i in pf:
            # if i.startswith(''):
            #     continue
            line=i.strip().split()
            c_centro_range=centro_region[line[0]]
            c_LTR_range=line[1:3]
            if has_intersection(c_LTR_range,c_centro_range):
                loc='C'
                C+=1
            else:
                OCIP+=1
                loc='OCIP'
                if 'Outperi' in line[3]:
                    loc='OC'
            outline='\t'.join(line)+'\t'+loc
            out.write(outline+'\n')
    print(C,OCIP)
if __name__ == "__main__":
    # pan gff txt
    # Chr01_SL187     RepeatMasker    dispersed_repeat        12841   13129   22.0    +       .       ID=2;Target "Motif:Helitron-N173_OS" 2265 2376
    # pan_gff_file_list=sys.argv[1]
    pericentro_LTR_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\merge_full_LTR.bed"
    # pan centro bed
    # Chr01_13-65 16749781 17724890
    pan_centro_bed_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed"
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_cluster\merge_full_LTR_LOCinfo.bed'
    main(pan_centro_bed_file,pericentro_LTR_bed_file,outfile)