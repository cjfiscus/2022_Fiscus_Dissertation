#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="calc"
#SBATCH -p koeniglab
#SBATCH --array=17-19

ARRAY="$SLURM_ARRAY_TASK_ID"

Rscript calc_r2.R ../results/"$ARRAY"mers.txt "$ARRAY" ../results/cor/"$ARRAY".txt
