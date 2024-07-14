# Script for non-redunant Te database build
  
## 05_repo_blast  
blast and parse repo sequence to calculate similarity  
  
## 06_gff_bed_new.sh  
terimal command used for comparisons between TE annotation based on repbase and merged TE database, to get TE entries exist over 10 and accumulate larger than 7.5k  
  
## 07_remove_dup.py  
remove redunant TE entries  in merged TE database (repbase + denovo)

## 08_repeatmask_rerun with 04 script  
batch repeatmask
  
## 09_full_stat    
script for full length and solo LTR identified  
