'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-31 21:39:39
LastEditors: Ne0tea
LastEditTime: 2024-02-01 13:24:18
'''
import sys
import os
def extract_sequences(fasta_file, bed_file,range):
    tmp_bed=open('tmp.bed','w')
    with open(bed_file, "r") as bed:
        for line in bed:
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            start_pos = max(0, start - int(range))
            end_pos = start
            record_id = f"{fields[3]}_forward"
            line='\t'.join([chrom,str(start_pos),str(end_pos),record_id,fields[4]])
            tmp_bed.write(line+'\n')

            start_pos = end
            end_pos =  end + int(range)
            record_id = f"{fields[3]}_backword"
            line='\t'.join([chrom,str(start_pos),str(end_pos),record_id,fields[4]])
            tmp_bed.write(line+'\n')
    cmd='bedtools getfasta -fi '+fasta_file+' -fo SZ-22_cluster_extract_fb'+str(range)+'.fasta -bed tmp.bed'
    os.system(cmd)
if __name__ == "__main__":
    #python 22_2_get_bed_fbseq.py .fasta .bed extract.fasta 500
    # fasta_file_path = "path/to/your/input.fasta"
    fasta_file_path=sys.argv[1]
    # bed_file_path = "path/to/your/input.bed"
    bed_file_path=sys.argv[2]
    # output_file_path = "path/to/your/output.fasta"
    range=sys.argv[3]
    extract_sequences(fasta_file_path, bed_file_path,range)
