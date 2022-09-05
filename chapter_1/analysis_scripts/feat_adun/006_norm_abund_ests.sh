#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --job-name="norm"
#SBATCH -p batch

# script to normalize K-mer based abundance using buscos

# assign filenames to variables
ABUND_IN="../../results/feat_abund/A_thal_kmer_abund.txt"
BUSCO_IN="../../results/feat_abund/A_thal_kmer_abund_busco.txt"
ABUND_NORM="../../results/feat_abund/A_thal_kmer_abund_norm.txt"
BUSCO_NORM="../../results/feat_abund/A_thal_kmer_iabund_busco_norm.txt"

# normalize
Rscript 006a_norm_seq_abund.R $ABUND_IN $BUSCO_IN $ABUND_NORM $BUSCO_NORM 
