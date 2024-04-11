'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-04-11 16:10:33
'''
import sys
import os
import logging
import copy
from itertools import chain
def find_longest_sublist(three_dim_list:dict):
    """
    寻找包含三维列表中的最长子列表（最内层），并返回包含该子列表的完整“列表组”。
    参数:
        three_dim_list (list): 三维列表，其中每个元素都是一个二维列表，每个二维列表的每个元素又是一个列表。
    返回:
        tuple: 一个包含最长子列表（最内层）和其所在“列表组”的元组。格式为 (longest_sublist, [outer_list_index, inner_list_index])。
                如果输入列表为空或没有子列表，则返回 (None, None)。
    """
    dic_list1=[x[0] for x in three_dim_list.values()]
    dic_list2=[x[1] for x in three_dim_list.values()]
    if 1:
        flat_list1 = list(chain.from_iterable(dic_list1))
        sorted_list1 = sorted(flat_list1, key=lambda x: int(x))
        sorted_list1 = list(set(sorted_list1))
        flat_list2 = list(chain.from_iterable(dic_list2))
        sorted_list2 = sorted(flat_list2, key=lambda x: int(x))
        sorted_list2 = list(set(sorted_list2))
        return sorted_list1,sorted_list2
    if 0:
        dic_listsum=[x[0]+x[1] for x in three_dim_list.values()]
        longest_sublist_index = max(((i, len(sublist)) for i, sublist in enumerate(dic_listsum)), key=lambda x: x[1])[0]
        # print(dic_list1[longest_sublist_index])
        return dic_list1[longest_sublist_index],dic_list2[longest_sublist_index]
    
def calculate_synteny_region_length(sp_pair,identity_dic,min_synteny_windowN=1):
    """
    计算并输出染色体共线性区域的长度。
    参数:
    sgbed (str): Stained glass bed 文件路径，用于分析共线性区域。
    min_synteny_windowN (int, 可选): 最小共线性窗口数，默认为5。
    返回:
    无显式返回值，但会输出共线性区域的信息到日志。
    """
    sp1=sp_pair[0]
    sp2=sp_pair[1]
    continuous_regions_dic={}
    dup_region_number=0
    dup_pair=[]
    # print(identity_dic)
    for i in identity_dic:
        continuous_regions=[]
        index=copy.deepcopy(i)
        for x in identity_dic[index]:
            continuous_regions2=[]
            index2=copy.deepcopy(x)
            c_pair1,c_pair2=[index,index2],[index2,index]
            if c_pair1 in dup_pair or c_pair2 in dup_pair:
                continue
            while index in identity_dic and index2 in identity_dic[index]:
                    # len(list(set(continuous_regions) & set(continuous_regions2))) <= 5:
                continuous_regions.append(index)
                continuous_regions2.append(index2)
                dup_pair.append([index,index2])
                index=str(int(index)+2000)
                index2=str(int(index2)+2000)
            '''same synteny region in different chr over 5 winodws was delete'''
            if len(continuous_regions) > min_synteny_windowN and len(continuous_regions2) > min_synteny_windowN:
                dup_region_number+=1
                continuous_regions_dic['Dup_number'+str(dup_region_number)]=[continuous_regions,continuous_regions2]
    '''output synteny region record'''
    # for i in continuous_regions_dic:
    #     loc1=continuous_regions_dic[i][0]
    #     loc2=continuous_regions_dic[i][1]
    #     for index_new in range(0,len(loc1)):
    #         num1,num2=loc1[index_new],loc2[index_new]
    #         loc=str(int(num1)-2000)+'_'+num1+'_'+str(int(num2)-2000)+'_'+num2
    #         # print(sg_dic[loc])
    #         # outfile.write(sg_dic[loc])
    if not continuous_regions_dic:
        return False,sp1,0,sp2,0
    LSR_chain_1,LSR_chain_2=find_longest_sublist(continuous_regions_dic)
    synteny_sp1,synteny_sp2=(int(len(LSR_chain_1))-1)*2000,(int(len(LSR_chain_2))-1)*2000

    logging.info(sp1+' Vs '+sp2+'\t'+','.join(LSR_chain_1)+'|'+','.join(LSR_chain_2))

    return True,sp1,synteny_sp1,sp2,synteny_sp2

if __name__ == "__main__":
    # SG_iden_bed=sys.argv[1]#Absulate path
    SG_iden_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_pw_SG\Chr01_YZ-2_CX3.data'
    centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\centromere_region_70m_used.bed'
    min_identity=95
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        filename='log_all_70_synteny_record',
                        filemode='w')
    centro_dic={}
    with open(centro_file,'r') as cf:
        for i in cf:
            line=i.strip().split()
            centro_dic[line[0]+'_centro']=int(line[2])-int(line[1])
    '''
    read stained glass bed file
    sg_dic* store bed file line as key
    identity_dic* store high identity region as chain
    '''
    try:
        last_chr1,last_chr2='',''
        with open(SG_iden_bed, 'r') as mf:
            sg_dic={}
            identity_dic={}
            synteny_len_dic={}
            for i in mf:
                line = i.strip().split()
                chr=line[0].split('_')[0]
                if i.startswith('#') or line[0] == line[3]:
                    continue
                if (last_chr1!=line[0] or last_chr2!=line[3]) and sg_dic and identity_dic:
                    condition,last_chr1,synteny_len1,last_chr2,synteny_len2=calculate_synteny_region_length((last_chr1,last_chr2),identity_dic)
                    if condition:
                        if last_chr1 in synteny_len_dic:
                            if last_chr2 in synteny_len_dic[last_chr1]:
                                synteny_len_dic[last_chr1][last_chr2]=synteny_len_dic[last_chr1][last_chr2]+synteny_len1/centro_dic[last_chr1]
                            else:
                                synteny_len_dic[last_chr1][last_chr2]=synteny_len1/centro_dic[last_chr1]
                        else:
                            synteny_len_dic[last_chr1]={last_chr2:synteny_len1/centro_dic[last_chr1]}
                    else:
                        pass
                    sg_dic={}
                    identity_dic={}
                if float(line[6]) > min_identity:
                    if line[2] in identity_dic:
                        identity_dic[line[2]].append(line[5])
                    else:
                        identity_dic[line[2]] = [line[5]]
                loc = '_'.join(line[1:3]) + '_' + '_'.join(line[4:6])
                sg_dic[loc] = i
                last_chr1,last_chr2=line[0],line[3]
            # print(identity_dic)
            condition,last_chr1,synteny_len1,last_chr2,synteny_len2=calculate_synteny_region_length((last_chr1,last_chr2),identity_dic)
            if last_chr1 in synteny_len_dic:
                if last_chr2 in synteny_len_dic[last_chr1]:
                    synteny_len_dic[last_chr1][last_chr2]=synteny_len_dic[last_chr1][last_chr2]+synteny_len1/centro_dic[last_chr1]
                else:
                    synteny_len_dic[last_chr1][last_chr2]=synteny_len1/centro_dic[last_chr1]
            else:
                synteny_len_dic[last_chr1]={last_chr2:synteny_len1/centro_dic[last_chr1]}
    except IOError as e:
        logging.error(f"File error: {e}")
        os.sys.exit(1)

    print(synteny_len_dic)
    print(last_chr1,synteny_len1,last_chr2,synteny_len2)
