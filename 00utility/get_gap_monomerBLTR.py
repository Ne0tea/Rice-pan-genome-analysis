'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-02-19 20:57:17
LastEditors: Ne0tea
LastEditTime: 2024-02-19 22:35:13
'''
def find_uncovered_ranges(range_list1:list, range_list2:list):
    start_point = range_list1[0][0]
    end_point = range_list1[-1][1]
    for i in range_list2[::-1]:
        if i[0]<start_point or i[1]>end_point:
            range_list2.remove(i)
    combined_ranges = sorted(range_list1 + range_list2, key=lambda x: x[0])
    uncovered_ranges = []
    # 检查是否有不连续的区域
    for i in range(len(combined_ranges) - 1):
        if combined_ranges[i][1] < combined_ranges[i+1][0]:
            uncovered_ranges.append((combined_ranges[i][1], combined_ranges[i+1][0]))

    return uncovered_ranges

def main(monomer,LTR):
    monomer_range={}
    with open(monomer,'r') as mf:
        for i in mf:
            line=i.strip().split()
            chr=line[0]
            start=int(line[1])
            end=int(line[2])
            if chr in monomer_range:
                last_end=monomer_range[chr][-1][1]
                if start==last_end+1:
                    monomer_range[chr][-1][1]=end
                else:
                    monomer_range[chr].append([start,end])
            else:
                monomer_range[chr]=[[start,end]]
    LTR_range={}
    with open(LTR,'r') as Lf:
        for i in Lf:
            line=i.strip().split()
            chr=line[0]
            start=int(line[1])
            end=int(line[2])
            if chr in LTR_range:
                LTR_range[chr].append([start,end])
            else:
                LTR_range[chr]=[[start,end]]
    uncover_region={}
    for i in monomer_range:
        chr_monomer=monomer_range[i]
        chr_LTR=LTR_range[i]
        uncover=find_uncovered_ranges(chr_monomer,chr_LTR)
        # print(uncover)
        uncover=[x for x in uncover if x[1]-x[0] > 2000]
        print(i,uncover)
        if not uncover:
            continue
        if i in uncover_region:
            uncover_region[i].append(uncover)
        else:
            uncover_region[i]=uncover
    # print(uncover_region)
    # for i in uncover_region:
    #     if i[1]-i[0] > 500:
    #         print(i)
if __name__ == "__main__":
    # monomer_file=sys.argv[1]#SZ-22_unit_centro.bed 
    monomer_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\monomer_anno.txt'
    # LTR_range_file=sys.argv[1]
    LTR_range_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\all_asm_peri_TE_70_replaced.gff3'
    main(monomer_file,LTR_range_file)