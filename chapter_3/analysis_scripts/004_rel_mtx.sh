#!/bin/bash
# calc relatedness mtx w/ gemma
# cjfiscus

# requires GEMMA 0.98.4

# define vars
GENO1=../data/cowpea_gwas
GENO2=../data/cowpea_envgwas
OUT=../results/gwas/

# calc centered relatedness matrix for gwas
#gemma -bfile "$GENO1" -gk 1 -outdir "$OUT" -o relmtx_gwas

# calc centered relatedness matrix for envgwas
gemma -bfile "$GENO2" -gk 1 -outdir "$OUT" -o relmtx_envgwas
