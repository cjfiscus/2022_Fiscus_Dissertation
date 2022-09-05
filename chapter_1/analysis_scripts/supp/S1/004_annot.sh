#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="annot"
#SBATCH -p koeniglab

# RepeatMasker 4.0.9

# this is TAIR10
GENOME=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

# annotate genome with library
RepeatMasker -gff -species arabidopsis -dir ../results/annotation "$GENOME"
