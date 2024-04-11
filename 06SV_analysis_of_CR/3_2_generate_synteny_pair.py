'''
Descripttion: gene_file:gene注释的bed文件
Author: Ne0tea
version:
Date: 2024-01-24 22:17:58
LastEditors: Ne0tea
LastEditTime: 2024-03-22 11:30:15
'''
import sys
import re
import logging
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def longest_common_substring(strings,accept_loss=0):
    if not strings:
        return []

    # 初始化最长公共子串
    longest_substring = []
    longest_substring_2nd = []
    longest_substring_3nd = []
    # 获取字符串的最短长度
    min_length = min(len(s) for s in strings)
    # 逐个字符比较
    for i in range(min_length):
        for j in range(min_length - i + 1):
            # 获取子串
            substring = strings[0][i:i+j]
            if sum(set(substring).issubset(set(s)) for s in strings) >= len(strings)-accept_loss:
                if len(substring) > len(longest_substring):
                    if not set(longest_substring).issubset(substring):
                        if not set(longest_substring_2nd).issubset(substring):
                            longest_substring_3nd=longest_substring_2nd
                        longest_substring_2nd=longest_substring
                    longest_substring = substring
                elif len(substring) > len(longest_substring_2nd) and not set(substring).issubset(longest_substring):
                    longest_substring_2nd = substring
                elif len(substring) > len(longest_substring_3nd) and not set(substring).issubset(longest_substring) \
                    and not set(substring).issubset(longest_substring_2nd):
                        longest_substring_3nd = substring
    return longest_substring,longest_substring_2nd,longest_substring_3nd

def has_intersection(rangeA, rangeB):
    # 分别提取两个范围的起始点和结束点
    start1, end1 = rangeA
    start2, end2 = rangeB
    c_has_inter=False
    if max(start1, start2) <= min(end1, end2):
        c_has_inter=True
    return c_has_inter

def main(gene_bed_file,accept_loss,phy_ord_file):
    phy_ord=[]
    with open(phy_ord_file,'r') as phyf:
        for i in phyf:
            line=i.strip().split()
            phy_ord.append(line[0])
    pattern=re.compile(r'LOC_Os(.*?)g')#LOC_Os10g
    gene_pattern=re.compile(r'ID=(.*?);Name')
    gene_bed_df=pd.read_table(gene_bed_file,names=['Chr','Start','End','ID','Strand','Location'])
    gene_bed_df[['Chr','Sp','_']]=gene_bed_df['Chr'].str.split('_', expand=True)
    for chr,group in gene_bed_df.groupby('Chr'):
        chr_df=pd.DataFrame()
        print(chr)
        # if chr!='Chr10':
        #     continue
        chr_centro_gene=[]
        Sp_list=[]
        grouped_group_df=group.groupby('Sp')
        for group_name in phy_ord:
            if group_name in grouped_group_df.groups:
                sp_group=grouped_group_df.get_group(group_name)
                sp=group_name
        # for sp,sp_group in group.groupby('Sp'):
                c_df=sp_group['ID'].str.split('|', expand=True)[1]
                c_df=c_df.str.replace(r'\.[^.]*\.[^.]*\.[^.]*$', '', regex=True)
                c_df=c_df.str.replace('fgenesh.mRNA.', '',regex=True)
                c_df=c_df.str.replace('LOC_Os', '',regex=True)
                c_df.name=sp

                last_upstream_index = sp_group[sp_group['Location'] == 'upstream'].index[-1]
                # print(sp_group)
                first_downstream_index = sp_group[sp_group['Location'] == 'downstream'].index[0]
                last_upstream_name=c_df[last_upstream_index]
                last_downstream_name=c_df[first_downstream_index]

                new_index_df = c_df.reset_index(drop=True)
                chr_df[sp]=new_index_df
                chr_centro_gene.append([last_upstream_name,last_downstream_name])
                Sp_list.append(sp)
        # if [id for id,x in chr_df.iteritems()] == Sp_list:
        #     print('Same')
        gene_list=[x.tolist() for id,x in chr_df.iteritems()]
        longest_substring,longest_substring_2nd,longest_substring_3nd=longest_common_substring(gene_list,accept_loss)
        print(longest_substring,longest_substring_2nd,longest_substring_3nd)

        plt.figure(figsize=(40, 5))
        draw_rect=[]
        texts=[]
        for id,i in enumerate(gene_list):
            left_site_x,left_site_y=23*id,0
            width=20
            height=0.22 * len(i)
            draw_rect.append([left_site_x,left_site_y,width,height])
            cur_texts=[]
            for x in i:
                if x not in longest_substring:
                    if x not in longest_substring_2nd:
                        if x not in longest_substring_3nd:
                            cur_texts.append((x,'none'))
                        else:
                            cur_texts.append((x,'orange'))
                    else:
                        cur_texts.append((x,'yellow'))
                else:
                    cur_texts.append((x,'red'))
            texts.append(cur_texts)
        # plt.gcf().subplots_adjust(left=0.01,top=0.91,bottom=0.09,right=0.98)
        for rect, text_list,centro_gene_list,Sp in zip(draw_rect, texts,chr_centro_gene,Sp_list):
            centro_start_gene=centro_gene_list[0]
            # if Sp=='CW14':
            #     print([x[0] for x in text_list])
            cur_sp_gene_list=[x[0] for x in text_list]
            centro_start_gene_index = len(cur_sp_gene_list) - 1 - cur_sp_gene_list[::-1].index(centro_start_gene)#centro区域起点基因为给定基因的最后一个匹配项
            centro_end_gene=centro_gene_list[1]
            centro_end_gene_index = cur_sp_gene_list.index(centro_end_gene)#centro区域终点基因为给定基因的第一个匹配项
            centro_gene_number=centro_end_gene_index-centro_start_gene_index

            # print(centro_start_gene,centro_end_gene)
            plt.gca().add_patch(patches.Rectangle((rect[0], rect[1]+0.22), rect[2], rect[3], color='lightgrey', alpha=0.3))
            plt.gca().add_patch(patches.Rectangle((rect[0]-1, (centro_start_gene_index+1)*0.22), 3, centro_gene_number*0.22, color='green', alpha=0.3))
            plt.text(rect[0]+10 , rect[1], Sp, fontsize=5,ha='center', va='center')
            for index,(text, bg_color) in enumerate(text_list):
                plt.text(rect[0]+10 , rect[1] + 0.22 * (index+1), text, fontsize=5,ha='center', va='center',
                        bbox=dict(facecolor=bg_color, alpha=0.5,edgecolor='none',pad=1))
        # plt.axis('equal')
        plt.axis('off')
        plt.xlim(0, 1625)
        plt.ylim(0, 10)
        # plt.show()
        # plt.invert_yaxis()
        plt.savefig('E:\\Bio_analysis\\Weedyrice\\pan_weedyrice\\gmap_CR\\'+chr+'_synteny.pdf', format='pdf', bbox_inches='tight')
        print('Draw plot finish')
        # chr_df.to_csv('E:\\Bio_analysis\\Weedyrice\\pan_weedyrice\\gmap_CR\\'+chr+'_gene_synteny.csv',index=False)

if __name__ == "__main__":
    # gene_bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\test.bed'
    gene_bed_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\gmap_CR\all_70_CR_20.filter.bed'
    phy_ord_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\70_phylogenic.txt'
    accept_loss=5
    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(gene_bed_file,accept_loss,phy_ord_file)
