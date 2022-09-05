#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="sim"
#SBATCH -p batch
#SBATCH --array=2-9

module load samtools/1.9

## define variables
COVERAGE=$(head -n $SLURM_ARRAY_TASK_ID ../data/coverage.txt | tail -n 1 | cut -f1)
READS=$(head -n $SLURM_ARRAY_TASK_ID ../data/coverage.txt | tail -n 1 | cut -f2)
REFERENCE=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

## simulate reads 
wgsim -N $READS -1 100 -2 100 -S9 "$REFERENCE" ../results/"$COVERAGE"_1.fq ../results/"$COVERAGE"_2.fq

## cleanup and make a single fq file for each sampling 
cat ../results/"$COVERAGE"_1.fq ../results/"$COVERAGE"_2.fq | gzip > ../results/"$COVERAGE"_sim.fq.gz
rm ../results/"$COVERAGE"_1.fq ../results/"$COVERAGE"_2.fq
