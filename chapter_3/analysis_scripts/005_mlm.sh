#!/bin/bash
# fit mlm with gemma
# cjfiscus
#SBATCH --job-name=gwas
#SBATCH -o ./std/%j.out
#SBATCH -e ./std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=2
#SBATCH --mem=16gb
#SBATCH -t 2-00:00:00
#SBATCH -p batch
#SBATCH --array=1-19

# requires GEMMA 0.98.5

# define variables 
GENO=../data/cowpea_envgwas # plink basename
KINSHIP=../results/gwas/relmtx_envgwas.cXX.txt # relatedness mtx
COL=$(($SLURM_ARRAY_TASK_ID)) # line in PHENOS and col (5 + N) in .fam 
PHENOS=../data/worldclim_vars.txt # list of phenotypes/climate vars
PHENO_NAME=$(head -n "$COL" "$PHENOS" | tail -n 1 | cut -f1) # parsed pheno/var name
OUT=../results/gwas

# association with mlm
echo "$PHENO_NAME"
#gemma -debug -bfile "$GENO" -n "$COL" -k $KINSHIP -lmm 4 -outdir "$OUT" -o $PHENO_NAME
gemma -debug -bfile "$GENO" -n "$COL" -c ../data/cov.txt -lm 4 -outdir "$OUT" -o $PHENO_NAME
