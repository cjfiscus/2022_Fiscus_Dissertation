#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100gb
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --time=1-00:15:00
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="filter"
#SBATCH -p koeniglab

# count tables
RAWCOUNTS="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers.txt"
FILTERED="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered.txt"

# apply initial filters
Rscript 004a_apply_filters.R "$RAWCOUNTS" "$FILTERED" 
