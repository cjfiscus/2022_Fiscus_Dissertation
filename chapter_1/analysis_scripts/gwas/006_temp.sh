#!/bin/bash

#SBATCH --job-name=post_gwas
#SBATCH -o std/%j.out
#SBATCH -e std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=4
#SBATCH --mem=16gb
#SBATCH -t 01:00:00
#SBATCH -p batch
#SBATCH --array=1-187

DIR=../../results/gwas/assoc2

NUM=$(echo "$SLURM_ARRAY_TASK_ID" + 0 | bc)
FILE=$(ls "$DIR"/*assoc.txt | head -n "$NUM" | tail -n1)
NAME=$(basename "$FILE" | sed 's/\.assoc\.txt//g')

# plot results 
Rscript 006a_plot_GWAS_gemma.R "$FILE" 
