#!/bin/bash
#SBATCH --job-name=corr_p
#SBATCH -o ./std/%j.out
#SBATCH -e ./std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=2
#SBATCH --mem=16gb
#SBATCH -t 1:00:00
#SBATCH -p koeniglab
#SBATCH --array=1-2384

WD=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/gwas/assoc
SCRIPT=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/scripts/gwas/008a_correct_p.R
NUM=$(echo "$SLURM_ARRAY_TASK_ID" + 2499 | bc )

cd $WD
FILE=$(ls *.assoc.txt | head -n "$NUM" | tail -n 1 | sed 's/*//g')

Rscript "$SCRIPT" "$FILE"
