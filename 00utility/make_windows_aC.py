'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-12-22 16:43:22
LastEditors: Ne0tea
LastEditTime: 2024-03-03 11:46:46
'''
import sys

def generate_regions(start, end, num_regions):
    chromosome_regions = []
    region_length = (end - start) // num_regions
    #print(region_length)
    step = region_length // 2
    for i in range(num_regions*2):
        region_start = max(start+step * i, 0)+1
        region_end = min(start+step * (i + 2), end)
        chromosome_regions.append((region_start, region_end))
    return chromosome_regions

def print_windows_info(chromosome, windows_list):
    for idx, window in enumerate(windows_list, start=1):
        start, end = window
        print(f"{chromosome} {idx} {start} {end}")

# 输出格式化的区域信息
def print_regions_info(chromosome, regions_list,startnumber):
    for idx, region in enumerate(regions_list, start=startnumber):
        start, end = region
        # print(f"{chromosome}_{idx}\t{start}\t{end}")
        print(f"{chromosome}\t{start}\t{end}")

def main(fai,centro_bed,peri_centro_bed):
    chr_length={}
    with open(fai,'r') as faif:
        for i in faif:
            line=i.strip().split()
            chr_length[line[0]]=int(line[1])
    sp_centro_dic={}
    with open(centro_bed,'r') as bf:
        for i in bf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            sp_centro_dic[line[0]]=[int(line[1]),int(line[2])]
    sp_pericentro_dic={}
    with open(peri_centro_bed,'r') as pcf:
        for i in pcf:
            line=i.strip().split()
            if line[0]=='Chr':
                continue
            sp_pericentro_dic[line[0]]=[int(line[1]),int(line[2])]
    for i in chr_length:
        chromosome_number = i  # 染色体号
        centromere_start = sp_centro_dic[i][0]  # 着丝粒起始位置
        centromere_end = sp_centro_dic[i][1]  # 着丝粒结束位置
        periC_start=sp_pericentro_dic[i][0]
        periC_end=sp_pericentro_dic[i][1]
        genome_size = chr_length[i]  # 基因组大小
        if centromere_start==centromere_end:
            continue
        num_regions = 20  # 每侧分割的区域数
        # print(periC_start,centromere_start,centromere_end,periC_end)
        # 计算着丝粒左侧的区域
        regions_L2periC = generate_regions(0, periC_start, num_regions)
        regions_periC2C = generate_regions(periC_start, centromere_start, num_regions)
        regions_C2C = generate_regions(centromere_start, centromere_end, num_regions)
        regions_C2periC = generate_regions(centromere_end, periC_end, num_regions)
        regions_periC2R = generate_regions(periC_end, genome_size, num_regions)

        # 0----------------PeriCentroStart----------------CentroStart----------------CentroEnd----------------PeriCentroEnd----------------End
        print_regions_info(chromosome_number, sorted(regions_L2periC, key=lambda x: x[0]), 1)
        print_regions_info(chromosome_number, sorted(regions_periC2C, key=lambda x: x[0]), num_regions+1)
        print_regions_info(chromosome_number, sorted(regions_C2C, key=lambda x: x[0]), num_regions*2+1)
        print_regions_info(chromosome_number, sorted(regions_C2periC, key=lambda x: x[0]), num_regions*3+1)
        print_regions_info(chromosome_number, sorted(regions_periC2R, key=lambda x: x[0]), num_regions*4+1)

    # print('Total bin number',num_regions*5*2)

if __name__ == "__main__":
    fai=sys.argv[1]#asm.fai made by samtools faidx
    # fai=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\CW03.fasta.fai'
    centro_bed=sys.argv[2]#bed file indicate centromere location chr1   Start   end
    # centro_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    peri_centro_bed=sys.argv[3]#bed file indicate centromere location chr1   Start   end
    # peri_centro_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_periregion_70m.bed'
    main(fai,centro_bed,peri_centro_bed)
