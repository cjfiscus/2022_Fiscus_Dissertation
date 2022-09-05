#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH -o ./std/%j.out
#SBATCH -e ./std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=2
#SBATCH --mem=16gb
#SBATCH -t 1-00:00:00
#SBATCH -p koeniglab
#SBATCH --array=1-42

# array 6-2499
# GEMMA 0.98.1

# define variables 
GENO=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/gwas/imputed3
KINSHIP=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/gwas/related_matrix.cXX.txt
PHENO_NAME=$(head -n "$SLURM_ARRAY_TASK_ID" /rhome/cfisc004/bigdata/projects/capsella_radseq/results/gwas/phenotypes.txt | tail -n 1)
COL=$(($SLURM_ARRAY_TASK_ID))
OUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/gwas/mlm

# association with mlm
echo "$PHENO_NAME"
gemma -debug -bfile "$GENO" -n "$COL" -k $KINSHIP -miss 0.20 -lmm 4 -outdir "$OUT" -o $PHENO_NAME
