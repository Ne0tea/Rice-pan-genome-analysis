'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-02-22 16:18:31
LastEditors: Ne0tea
LastEditTime: 2024-02-23 15:02:30
'''
import pandas as pd
import argparse

def get_range_meth(df_all,df_all_WGBS,start,end):
    df_freq = df_all[(df_all["start"] > start) & (df_all["end"] <= end)].reset_index()
    df_chr_WGBS = df_all_WGBS[(df_all_WGBS["pos"] > start) & (df_all_WGBS["pos"] <= end)].reset_index(drop=True)
    if df_freq.empty:
        prope = 0
    else:
        C_count = df_chr_WGBS.shape[0]
        prope = df_freq["probability"].sum() / C_count
    return prope

def get_centromere_region(prefix, chr, outfile,target_range,windows=10000):
    df_all = pd.read_table("/public2/labmember/xielj/T2Trice_centromere/06_Methylation/" + str(chr) + "_"  + str(prefix) + "_ccs.align.sorted.bam_to_cpg.combined.bed", sep='\t', header=None)
    df_all.columns = ["chr", "start", "end", "score", "haplotype", "coverage", "modified", "unmodified", "probability"]

    df_all_WGBS = pd.read_table("/public2/labmember/xielj/T2Trice_centromere/00_data/methylation_data/Bismark/" + str(prefix) + "_result/" + str(prefix) + "_report/" + str(prefix) + "." + str(prefix) + "_1_paired_bismark_bt2_pe.deduplicated.CX_report.txt.chr" + str(chr) + "_" + str(prefix) + ".CX_report.txt.gz", sep='\t', header=None)
    df_all_WGBS.columns = ["chr", "pos", "strand", "modified", "unmodified", "type", "nucl"]

    target_range_df=pd.read_table(target_range, sep='\t', header=None,usecols=[0,1,2,3])
    target_range_df.columns = ["chr", "start", "end","unit_name"]
    target_range_df = target_range_df.sort_values(by='start')
    o = open(outfile, 'a')
    ful_s=0
    ful_e=0+windows
    for index,row in target_range_df.iterrows():
        s=row['start']
        e=row['end']
        while ful_e < s:
            prope=get_range_meth(df_all,df_all_WGBS,ful_s,ful_e)
            o.write(str(chr) + "_" + str(prefix) + "\t" + str(ful_s) + "\t" + str(ful_e) + "\t" + str(prope) + "\tOUT_Target\n")
            ful_s,ful_e=ful_s+windows,ful_e+windows
        prope=get_range_meth(df_all,df_all_WGBS,ful_s,s)
        o.write(str(chr) + "_" + str(prefix) + "\t" + str(ful_s) + "\t" + str(s) + "\t" + str(prope) + "\tOUT_Target\n")
        prope=get_range_meth(df_all,df_all_WGBS,s,e)
        o.write(str(chr) + "_" + str(prefix) + "\t" + str(s) + "\t" + str(e) + "\t" + str(prope) + "\tTarget\n")
        ful_s,ful_e=e,e+windows
    o.close()
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='prefix', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='chr', dest='input2', required=True)
    parser.add_argument('-range', '--range', help='targetR', dest='range', required=True)
    parser.add_argument('-windows', '--windows', help='win', dest='win',default=10000)
    # parser.add_argument('-r2', '--range2', help='chr', dest='input2', required=True)
    parser.add_argument('-o', '--output', help='out', dest='output', required=True)
    args = parser.parse_args()
    get_centromere_region(args.input1, args.input2, args.output,args.range,args.win)

if __name__ == '__main__':
    main()
