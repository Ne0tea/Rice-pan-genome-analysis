'''
Descripttion: restore costume method
Author: Ne0tea
version: 
Date: 2024-04-02 15:08:01
LastEditors: Ne0tea
LastEditTime: 2024-04-02 17:11:05
'''
def has_intersection(rangeA, rangeB):
    '''
    判断两个范围是否存在交集且不完全重叠。

    参数:
    rangeA: 第一个范围，表示为一个包含起始点和结束点的元组。
    rangeB: 第二个范围，同样表示为一个包含起始点和结束点的元组。

    返回值:
    返回一个布尔值，如果两个范围存在交集且不完全重叠，则为True；否则为False。
    '''
    # 提取范围的起始点和结束点
    start1, end1 = rangeA
    start2, end2 = rangeB
    c_has_inter = False  # 初始化没有交集的标志

    # 判断两个范围是否存在交集且不完全重叠
    if max(start1, start2) <= min(end1, end2) and not (start1 < start2 and end1 > end2):
        c_has_inter = True

    return c_has_inter

def list_range_include(sort_range_list,target_range):
    """
    判断目标范围是否与给定排序后的范围列表中的某个范围完全重叠。
    :param sort_range_list: 排序后的范围列表，每个范围由起始和结束点组成。
    :param target_range: 目标范围，同样由起始和结束点组成。
    :return: 布尔值，如果目标范围与列表中某个范围完全重叠，则返回True；否则返回False。
    """
    cur_range_list=sort_range_list[:]
    cur_range_list.append(target_range)
    '''
    分别将ref列表根据范围起点排序，提取target在ref中左右两侧的范围，判断是否有完全重叠
    '''
    cur_range_list_sort=sorted(cur_range_list, key=lambda x: x[0])
    target_index=cur_range_list_sort.index(target_range)
    if target_index == 0:
        right_range=cur_range_list_sort[target_index+1]
        if right_range[0] <= target_range[0]+50 and right_range[1] >= target_range[1]-50:
            return True
    else:
        if target_index == len(cur_range_list_sort)-1:
            left_range=cur_range_list_sort[target_index-1]
            if left_range[0] <= target_range[0]+50 and left_range[1] >= target_range[1]-50:
                return True
        else:
            left_range=cur_range_list_sort[target_index-1]
            right_range=cur_range_list_sort[target_index+1]
            for i in [left_range,right_range]:
                if i[0] <= target_range[0]+50 and i[1] >= target_range[1]-50:
                    return True
    return False

def intersection_length(range1, range2):
    start1, end1 = range1
    start2, end2 = range2
    # 计算交集的起始点和结束点
    start = max(start1, start2)
    end = min(end1, end2)
    # 如果交集存在，则返回交集长度；否则返回0
    if start <= end:
        return end - start + 1
    else:
        return 0