#!/bin/bash

#SBATCH --job-name=xpclr
#SBATCH -o std/%j.out
#SBATCH -e std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=4
#SBATCH --mem=64gb
#SBATCH -t 1-00:00:00
#SBATCH -p koeniglab
#SBATCH --array=1-11

# xpclr
A=../data/breeding.txt
B=../data/wild.txt
I="../data/plink.vcf"
F="vcf"
O="../results/xpclr"
C_LST=../data/chr_lst.txt
C=$(head -n "$SLURM_ARRAY_TASK_ID" "$C_LST" | tail -n 1)

source activate xpclr

xpclr -O "$O"/"$C"_wildvbreeding_xpclr.txt -F "$F" -I "$I" -Sa "$A" -Sb "$B" -C "$C" --verbose 10 --minsnps 2
