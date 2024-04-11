'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-11 23:08:27
LastEditors: Ne0tea
LastEditTime: 2023-11-29 10:57:03
'''
import numpy as np
import sys

def calculate_standard(data):
    # 获取数据的属性
    attr1 = [data[d][0]+data[d][1]for d in data]  # 第一个属性
    # attr2 = [data[d][1] for d in data]  # 第二个属性

    # 计算平均值和标准差
    mean_attr1 = np.mean(attr1)
    std_attr1 = np.std(attr1)

    # mean_attr2 = np.mean(attr2)
    # std_attr2 = np.std(attr2)

    # 根据正态分布计算标准
    # 对于每个属性，标准为 (均值 - 1.96 * 标准差, 均值 + 1.96 * 标准差)
    std_attr1_lower = mean_attr1 - 1.96 * std_attr1
    std_attr1_upper = mean_attr1 + 1.96 * std_attr1

    # std_attr2_lower = mean_attr2 - 1.96 * std_attr2
    # std_attr2_upper = mean_attr2 + 1.96 * std_attr2

    # 返回标准
    # return [(std_attr1_lower, std_attr1_upper), (std_attr2_lower, std_attr2_upper)]
    return  (std_attr1_lower, std_attr1_upper)

def main(blast):
    self_blast_dic={}
    with open(blast,'r') as blastf:
        for i in blastf:
            line=i.strip().split()
            ratio=float(line[5]) if float(line[5]) < 1 else 1
            self_blast_dic[line[0]]=[float(line[2]),ratio*100]
    # print(self_blast_dic)
    createrial=calculate_standard(self_blast_dic)
    print(createrial)
if __name__ == "__main__":
    # blast_file=sys.argv[1]
    blast_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\oryrep.ref_blast.out_parse'
    main(blast_file)