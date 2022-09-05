#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="gc"
#SBATCH --time=1-00:15:00
#SBATCH -p batch
#SBATCH --array=2-1172

KMER_TABLE="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered_gc.txt"
OUT="../results/gc/parts"

Rscript 006a_corr_count_gc.R "$KMER_TABLE" "$OUT"/"$SLURM_ARRAY_TASK_ID".txt "$SLURM_ARRAY_TASK_ID"
