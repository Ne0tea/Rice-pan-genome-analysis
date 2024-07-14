'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-02-13 21:18:25
'''
import sys
import logging

def main(monomer_file,ID_loc_file):
    ID_monomer_dic={}
    ID_loc_dic={} 
    with open(ID_loc_file,'r') as IDf:
        for i in IDf:
            line=i.strip().split()
            ID_loc_dic[line[3]+'_fw']=[int(line[1]),int(line[2])]
    with open(monomer_file,'r') as monof:
        for i in monof:
            line=i.strip().split()
            ID_monomer_dic[line[0]]=line[1]
    _dict = {}
    for key, value in ID_monomer_dic.items():
        if not value:
            continue
        if value not in _dict.keys():
            _dict[value] = []
        _dict[value].append(key)
    
    for key, value in sorted(_dict.items(), key=lambda x: len(x[1]),reverse=True):
        if len(value) > 1:
            _str = ",".join([str(x) for x in value])
            # print("Same key in dic {key} is {values}, location is {loc}".format(key=key, values=_str,loc=[ID_loc_dic[x] for x in value]))
            print(key,len(value),_str,[ID_loc_dic[x] for x in value],sep='\t')
if __name__ == "__main__":
    monomer_file=sys.argv[1]
    # monomer_file=r'D:\study resource\Digitaria\Script\SZ-22_ID_monomer.txt'
    ID_loc_file=sys.argv[2]
    # ID_loc_file=r'D:\study resource\tmp.bed'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(monomer_file,ID_loc_file)
