#!/bin/bash

#SBATCH --job-name=post_gwas
#SBATCH -o std/%j.out
#SBATCH -e std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=4
#SBATCH --mem=16gb
#SBATCH -t 01:00:00
#SBATCH -p intel
#SBATCH --array=1-1373

NUM=$(echo "$SLURM_ARRAY_TASK_ID" + 0 | bc)

DIR=../../results/gwas/assoc2
FILE=$(ls "$DIR"/*assoc.txt | head -n "$NUM" | tail -n1)
NAME=$(basename "$FILE" | sed 's/\.assoc\.txt//g')

# determine significance threshold
THRESHOLD=0.00001

# print sig snps 
cat "$FILE" | awk -v THRESHOLD="$THRESHOLD" -v NAME="$NAME" '$14 < THRESHOLD {print NAME,"\t",$0}' "$FILE1" \
	| sed 's/.assoc.txt//g' > "$DIR"/"$NAME"_snps_lowp.txt
