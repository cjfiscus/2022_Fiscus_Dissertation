#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH -o ./std/%j.out
#SBATCH -e ./std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=2
#SBATCH --mem=16gb
#SBATCH -t 1-00:00:00
#SBATCH -p intel
#SBATCH --array=2382-2384

# GEMMA 0.98.1

# define variables 
GENO=../../results/gwas/genos_filtered
KINSHIP=../../results/gwas/related_matrix.cXX.txt
NUM=$(echo "$SLURM_ARRAY_TASK_ID" + 2499 | bc)
PHENO_NAME=$(head -n "$NUM" ../../results/gwas/pheno_lst.txt | tail -n 1)
COL=$(($NUM))
OUT=../../results/gwas/assoc

# association with mlm
echo "$PHENO_NAME"
gemma -debug -bfile "$GENO" -n "$COL" -k $KINSHIP -lmm 4 -outdir "$OUT" -o $PHENO_NAME
