'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-14 20:17:28
LastEditors: Ne0tea
LastEditTime: 2024-01-20 16:55:58
'''
import matplotlib.pyplot as plt
import numpy as np
import os

identity_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\SZ-22_identity.line'
OC_identity_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\SZ-22\OUTSZ-22_identity.line'
filedir=os.path.dirname(identity_file)
filename=os.path.basename(identity_file)
iden_list=[]
with open(identity_file,'r') as ifile:
    iden_list=ifile.readline().strip().split()[1:]
iden_list = [float(value) for value in iden_list]

OCiden_list=[]
with open(OC_identity_file,'r') as ifile:
    OCiden_list=ifile.readline().strip().split()[1:]
OCiden_list = [float(value) for value in OCiden_list]

iden_list = np.asarray(iden_list, dtype='float64')
bin_edges = np.linspace(min(iden_list), max(iden_list), num=100)
# 绘制直方图
plt.hist(OCiden_list, bins=bin_edges, alpha=0.5, label='Chromosome arm', color='#283618', edgecolor='none')
plt.hist(iden_list, bins=bin_edges, alpha=0.5, label='Centric',color='#bc6c25', edgecolor='none')# 调整 bins 的数量以更改直方图的精细程度

plt.title('Histogram')
plt.xlabel('Values')
plt.ylabel('Frequency')

x_ticks = np.linspace(min(min(iden_list), min(OCiden_list)), max(max(iden_list), max(OCiden_list)), 9)  # 生成9个分割符
plt.xticks(x_ticks[::2])
# plt.savefig(filedir+'/'+filename+'.jpg')
# 显示图形
plt.show()
