#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="samp"
#SBATCH -p koeniglab

# sample raw counts
IN=/rhome/cfisc004/bigdata/projects/K-mers_Arabidopsis/results/kmer_counts/A_thaliana_12mers.txt.gz
OUT=~/bigdata/projects/SX/results/counts/1741_Raw.txt.gz
zcat "$IN" | cut -f1,69 | tail -n+2 | gzip > "$OUT"

# sample normalized counts 
#IN=/rhome/cfisc004/bigdata/projects/K-mers_Arabidopsis/results/kmer_counts/A_thaliana_12mers_anno_GC_TMM_filtered.txt
#OUT=~/bigdata/projects/SX/results/counts/1741_Norm.txt.gz
#cut -f2,50 "$IN" | tail -n+2 | gzip > "$OUT"
