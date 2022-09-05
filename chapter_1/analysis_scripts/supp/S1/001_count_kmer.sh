#!/bin/bash -l 

#SBATCH --job-name=mers
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=2
#SBATCH --mem=16gb
#SBATCH -t 02:00:00
#SBATCH -p koeniglab
#SBATCH --array=12

module load jellyfish/2.2.9

# parameters set for Col-0
REFERENCE=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.chromosomes_only.fa

## count kmers in ref genome
jellyfish count -C -m "$SLURM_ARRAY_TASK_ID" -s 500M -t 2 -o ../results/TAIR10_"$SLURM_ARRAY_TASK_ID".jf "$REFERENCE"
jellyfish dump -tc ../results/TAIR10_"$SLURM_ARRAY_TASK_ID".jf | gzip > ../results/TAIR10_"$SLURM_ARRAY_TASK_ID".txt.gz

# cleanup
rm ../results/*.jf
