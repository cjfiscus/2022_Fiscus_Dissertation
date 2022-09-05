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
Rscript X_plot_GWAS_gemma_loop2.R plot_lst4.txt
