'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-04-05 23:01:32
'''
import sys
import logging
import copy

def rm_dup(lists):
    """
    去除列表中的重复元素，但是保留最长的重复元素列表。
    参数:
    lists -- 包含多个子列表的列表
    返回值:
    unique_lists -- 去除重复元素后的列表，其中如果存在重复元素，则保留最长的重复元素列表
    """
    unique_lists = []
    for index1 in range(0,len(lists)):
        homo_list=[]
        homo_list.append(lists[index1])
        for index2 in range(0,len(lists)):
            if len(list(set(lists[index1]) & set(lists[index2]))) != 0:
                homo_list.append(lists[index2])
        if len(homo_list)==0:
            unique_lists.append(lists[index1])
        longest_list = max(homo_list, key=len)
        if longest_list == lists[index1]:
            unique_lists.append(longest_list)
    return unique_lists

def main(sgbed):
    outfile=open(SG_iden_bed.split('.')[0]+'_Continue.bed','w')
    sg_dic={}
    identity_dic={}
    with open(sgbed,'r') as mf:
        for i in mf:
            if i.startswith('#'):
                continue
            line=i.strip().split()
            if float(line[6]) > 93:
                if line[2] in identity_dic:
                    identity_dic[line[2]].append(line[5])
                else:
                    identity_dic[line[2]]=[line[5]]
            loc='_'.join(line[1:3])+'_'+'_'.join(line[4:6])
            sg_dic[loc]=i
    # print(identity_dic)
    continuous_regions_dic={}
    dup_region_number=0
    dup_pair=[]
    for i in identity_dic:
        continuous_regions=[]
        index=copy.deepcopy(i)
        for x in identity_dic[index]:
            continuous_regions2=[]
            index2=copy.deepcopy(x)
            c_pair1,c_pair2=[index,index2],[index2,index]
            if c_pair1 in dup_pair or c_pair2 in dup_pair:
                continue
            while index2 in identity_dic and index in identity_dic and index2 in identity_dic[index] and \
                    len(list(set(continuous_regions) & set(continuous_regions2))) <= 5:
                continuous_regions.append(index)
                continuous_regions2.append(index2)
                dup_pair.append([index,index2])
                index=str(int(index)+2000)
                index2=str(int(index2)+2000)

            if len(continuous_regions) > 15 and len(continuous_regions2) > 15:
                dup_region_number+=1
                continuous_regions_dic['Dup_number'+str(dup_region_number)]=[continuous_regions,continuous_regions2]

    # print(continuous_regions_dic)

    for i in continuous_regions_dic:
        loc1=continuous_regions_dic[i][0]
        loc2=continuous_regions_dic[i][1]
        for index_new in range(0,len(loc1)):
            num1,num2=loc1[index_new],loc2[index_new]
            loc=str(int(num1)-2000)+'_'+num1+'_'+str(int(num2)-2000)+'_'+num2
            print(sg_dic[loc])
    #         outfile.write(sg_dic[loc])
    # outfile.close()
    for i in continuous_regions_dic:
        logging.info(line[0]+'\t'+','.join(continuous_regions_dic[i][0]))

if __name__ == "__main__":
    # SG_iden_bed=sys.argv[1]#Absulate path
    SG_iden_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_pw_SG\Chr01_NJ11_centro.2000.5000.bed'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(SG_iden_bed)
