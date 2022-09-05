#!/bin/bash
# make datasets for gwas and envgwas
# cjfiscus
# 2021-11-09
# updated 2022-03-28

# requires plink 1.9
plink --tfile ../data/cowpea --allow-extra-chr --out ../data/cowpea_gwas

# copies of binary files for egwas analyses
cp ../data/cowpea_gwas.bed ../data/cowpea_envgwas.bed
cp ../data/cowpea_gwas.bim ../data/cowpea_envgwas.bim
