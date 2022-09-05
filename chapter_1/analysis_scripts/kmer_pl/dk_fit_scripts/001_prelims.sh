cd /rhome/dkoenig/bigdata/SCRATCH/MODEL_FIT_SEQ_BY_REMOVE/DATA/OUTPUT
mkdir SPLIT_IT
cd SPLIT_IT
mkdir ../CORRECTED
split -l 10000 /rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered_gc_qq_final.txt
head -n1 xaa > ../CORRECTED/header.txt
grep -v mer xaa >temp
mv temp xaa
