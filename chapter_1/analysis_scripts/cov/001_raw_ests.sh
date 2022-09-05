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
DIR=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/genome
INPUT=$(ls "$DIR"/*.per-base.bed.gz | head -n "$SLURM_ARRAY_TASK_ID" | tail -n1)
NAME=$(basename "$INPUT" | cut -f1 -d".")
REPEATBED=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/data/TAIR10_repeats_sorted.bed
BUSCOBED=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/data/A_thal_complete_buscos.bed
OUT1=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/ests/"$NAME"_repeats.txt
OUT2=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/ests/"$NAME"_busco.txt

## total coverage across all repeats
bedtools map -a "$REPEATBED" -b "$INPUT" -c 4 -o median > "$OUT1"

## median coverage per busco annotation
bedtools map -a "$BUSCOBED" -b "$INPUT" -c 4 -o median > "$OUT2"
