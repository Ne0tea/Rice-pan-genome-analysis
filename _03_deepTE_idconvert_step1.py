'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-09 11:38:23
LastEditors: Ne0tea
LastEditTime: 2023-11-11 23:09:51
'''
import sys
Class_table_file=sys.argv[1]
id_file=sys.argv[2]
out_file=open('deepTE_idconverted.id','w')
# Class_table_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\Class_convert_table.txt"
# id_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\opt_DeepTE.txt'
# id_file=r"E:\Bio_analysis\Weedyrice\pan_weedyrice\oryrep.id"

type_list=[]
class_dic={}
with open(Class_table_file,'r') as classf:
    for i in classf:
        line=i.strip().split('\t')
        class_dic[line[0]]=line[1]

with open(id_file,'r') as seqf:
    if id_file.endswith('oryrep.id'):
        for i in seqf:
            line=i.strip().split('|')
            order=class_dic[line[1]]
            print(line[0]+'#'+order)
    else:
        for i in seqf:
            line=i.strip().split()
            order=class_dic[line[1]]
            repeatmodeler_name=line[0].split('#')[0]
            repeatmodeler_class=line[0].split('#')[1]
            if order=='Unknown':
                out_file.write(repeatmodeler_name+'#'+repeatmodeler_class+'\n')
            elif repeatmodeler_class =='Unknown' and order!='Unknown':
                out_file.write(repeatmodeler_name+'#'+order+'\n')
            elif repeatmodeler_class =='LTR/Unknown' and order!='Unknown':
                out_file.write(repeatmodeler_name+'#'+order+'\n')
            elif repeatmodeler_class != 'Unknown' and repeatmodeler_class !='LTR/Unknown':
                out_file.write(repeatmodeler_name+'#'+repeatmodeler_class+'\n')
            else:
                print(repeatmodeler_name,'wasn\'t fully pasesed')
print(id_file.split('/')[-1],'was converted!')