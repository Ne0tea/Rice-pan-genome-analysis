'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-30 15:07:26
LastEditors: Ne0tea
LastEditTime: 2024-01-30 15:41:25
'''

def main(pericentro_region_file,unit_bed,outfile):
    pass
    peri_dic={}
    with open(pericentro_region_file,'r') as prf:
        for i in prf:
            line=i.strip().split()
            peri_dic[line[0]]=[int(line[1]),int(line[2])]
    of=open(outfile,'w')
    number=0
    with open(unit_bed,'r') as ubf:
        for i in ubf:
            number+=1
            line=i.strip().split()
            c_pericentro_region_s=peri_dic[line[0]][0]
            c_pericentro_region_e=peri_dic[line[0]][1]
            unit_s=int(line[1])
            unit_e=int(line[2])
            old_name=line[3]
            #new name Chr08_NJ11_SZ-22_LTR_centro_Cluster237_C_2
            #new name Chr08_NJ11_SZ-22_LTR_centro_intact_4_1
            if unit_e < c_pericentro_region_s or unit_s > c_pericentro_region_e:
                location='Outperi'
            else:
                location='Inperi'
            seq_name='_'.join(line[3].split('_')[:3])+'_intact_'+str(number)
            outline='\t'.join(line[:3])+'\t'+seq_name+'\t'+'\t'.join(line[4:6])+'\t'+location
            of.write(outline+'\n')
    of.close()
if __name__ == "__main__":
    pericentro_region_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_periregion_70m.bed"
    unit_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\Full_SZ-22_intact.bed'
    outfile=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\Full_SZ-22_intact_peri.bed'
    main(pericentro_region_file,unit_bed,outfile)