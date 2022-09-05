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
#SBATCH --array=1-540

# array 6-2499
# GEMMA 0.98.1

# define variables 
GENO=../../results/gwas/genos_filtered4
KINSHIP=../../results/gwas/related_matrix.cXX.txt
PHENO_NAME=$(head -n "$SLURM_ARRAY_TASK_ID" ../../results/gwas/phenotypes_cov.txt | tail -n 1)
COL=$(($SLURM_ARRAY_TASK_ID))
OUT=../../results/gwas/assoc3

# association with mlm
echo "$PHENO_NAME"
gemma -debug -bfile "$GENO" -n "$COL" -k $KINSHIP -lmm 4 -outdir "$OUT" -o $PHENO_NAME
