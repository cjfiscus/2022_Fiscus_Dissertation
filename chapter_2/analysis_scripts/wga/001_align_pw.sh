#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=64G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="aln"
#SBATCH -p koeniglab
#SBATCH --array=3-5

## dependency mummer 4.0rc1
MAP=../../data/alignments_pw.txt
OUT=../../results/alignments/pairwise
NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$MAP" | tail -n1 | cut -f1)
GEN1=$(head -n "$SLURM_ARRAY_TASK_ID" "$MAP" | tail -n1 | cut -f4)
GEN2=$(head -n "$SLURM_ARRAY_TASK_ID" "$MAP" | tail -n1 | cut -f5)

# align with mummer (produce both delta and sam)
echo "aligning "$GEN1" to "$GEN2""
nucmer -t 4 --sam-long="$OUT"/"$NAME" "$GEN1" "$GEN2"
nucmer -t 4 -p "$OUT"/"$NAME" "$GEN1" "$GEN2"
