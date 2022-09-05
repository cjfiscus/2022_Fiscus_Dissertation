#!/bin/bash

#SBATCH --job-name=ld
#SBATCH -o std/%j.out
#SBATCH -e std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=4
#SBATCH --mem=32gb
#SBATCH -t 3-00:00:00
#SBATCH -p koeniglab

module load plink/1.90b6.25

IN=/rhome/cfisc004/bigdata/projects/cowpea_gwas/data/cowpea_gwas

plink --memory 32000 --threads 4 --bfile "$IN" --allow-extra-chr --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 10000 --out ld

