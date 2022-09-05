#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="cov"
#SBATCH -p koeniglab
#SBATCH --array=1-1269

## array 1-1269

# software dependencies
## bedtools 2.27.1; 
conda deactivate

# SET VARIABLES
DIR=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/ests
INPUT=$(ls "$DIR"/*_repeats.txt | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1)
NAME=$(basename "$INPUT" | sed 's/_repeats\.txt//g')

Rscript 002b_norm_ests.R "$NAME"
