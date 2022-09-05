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
#SBATCH --array=1-16

# xpclr
A=../data/cluster1.txt
B=../data/cluster2.txt
I="/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops3/populations.snps.vcf"
F="vcf"
O="../results/xpclr"
C_LST=../data/chr_lst.txt
C=$(head -n "$SLURM_ARRAY_TASK_ID" "$C_LST" | tail -n 1)

source activate xpclr

xpclr -O "$O"/"$C"_xpclr.txt -F "$F" -I "$I" -Sa "$A" -Sb "$B" -C "$C" --verbose 10 --minsnps 2 --size 100000
