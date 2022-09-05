#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="sm"
#SBATCH -p koeniglab
#SBATCH --array=10

# RepeatMasker 4.1.1
# using RMBlast 2.10.0
# using Dfam 3.2
# RepBase RepeatMasker Edition 20181026

module unload miniconda2
module unload python perl
module load miniconda3
conda activate my-perl-env

# this is TAIR10
LST=../../data/genomes.txt
GENOME=$(head -n "$SLURM_ARRAY_TASK_ID" "$LST" | tail -n1)

# annotate genome with library
RepeatMasker -gff -xsmall -dir ../../results/annotation "$GENOME"
