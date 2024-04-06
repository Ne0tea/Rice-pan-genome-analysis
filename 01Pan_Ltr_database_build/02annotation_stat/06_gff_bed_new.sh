echo "Starting job at:"
date -d +20minute +'%Y-%m-%d %H:%M:%S'
###job dir

###variable
#query=""
repbase_mask_dir="/public/home/huangyj/pan_rice/${1}_repbase"
repdeep_mask_dir="/public/home/huangyj/pan_rice/${1}_repdeep"
rep_base="/public/home/huangyj/pan_rice/oryrep.ref"
###job command
#source activate /public/home/huangyj/anaconda3/envs/r_p

mkdir -p gff_bed_new && cd gff_bed_new
awk 'OFS="\t"{print$1,$4,$5,$10}' ${repbase_mask_dir}/${1}.fasta.out.gff > ${1}_repbase.bed
awk 'OFS="\t"{print$1,$4,$5,$10}' ${repdeep_mask_dir}/${1}.fasta.out.gff > ${1}_repdeep.bed
bedtools subtract -A -a ${1}_repdeep.bed -b ${1}_repbase.bed > ${1}_repdeep_rep_sub.bed
python /public/home/huangyj/pan_rice/make_pan/blast_out/cal_len_per_seq_from_bed.py ${1}_repdeep_rep_sub.bed 10 7500
ls /public/home/huangyj/pan_rice/make_pan/blast_out/*${1}*_parse  > ${1}_parsed.list
python /public/home/huangyj/pan_rice/make_pan/blast_out/get_new_seq.py ${1}_parsed.list ${1}_repdeep_rep_sub.bedcontribution.txt
cd ../
#conda deactivate
#conda deactivate
echo "Finished at:"
date -d +20minute +'%Y-%m-%d %H:%M:%S'

