#!/bin/bash

#SBATCH --job-name=plot
#SBATCH -o std/%j.out
#SBATCH -e std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=4
#SBATCH --mem=64gb
#SBATCH -t 1-00:00:00
#SBATCH -p koeniglab

# plot results
ls /rhome/cfisc004/bigdata/projects/capsella_radseq/results/gwas/mlm/*.assoc.txt > plot_lst.txt

Rscript 016a_plot_GWAS_gemma_loop.R plot_lst.txt capsella_gwas.pdf
