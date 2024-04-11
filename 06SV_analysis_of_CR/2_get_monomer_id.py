'''
Descripttion: file1:monomer注释文件;file2:gene注释文件;筛选在gene注释末端的monomer序列
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-03-18 14:39:15
'''
import sys
import os
import pandas as pd
import logging

def main(monomer_bed,wkdir):
    data = pd.read_table(monomer_bed,names=['Chr','Start','End','name','Strand','group'])
    type_counts = data.iloc[:, 5].value_counts()

    # 提取出现最多的五种类型
    top_five_types = type_counts.head(5).index.tolist()
    # 输出包含这五种类型的行
    rows_with_top_types = data[data.iloc[:, 5].isin(top_five_types)]
    os.makedirs(wkdir, exist_ok=True)
    for type_name in top_five_types:
        rows_with_type = data[data.iloc[:, 5] == type_name]
        file_name = os.path.join(wkdir, f"{type_name}_idname.bed")
        rows_with_type.to_csv(file_name, index=False,sep='\t',header=False)
        print(f"Saved {len(rows_with_type)} rows with type '{type_name}' to '{file_name}'")

if __name__ == "__main__":
    wkdir=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR'
    monomer_bed=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\Chr08_monomer_besideG.bed'
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(monomer_bed,wkdir)
