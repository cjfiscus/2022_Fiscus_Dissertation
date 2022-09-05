#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="samp"
#SBATCH -p koeniglab
#SBATCH --array=9

module load seqtk/281537b

COVERAGE=$(head -n $SLURM_ARRAY_TASK_ID ../data/coverage.txt | tail -n 1 | cut -f1)
READS=$(head -n $SLURM_ARRAY_TASK_ID ../data/coverage.txt | tail -n 1 | cut -f2)

# sample reads
seqtk sample -s9 ../data/6909_1.fastq.gz $READS > ../results/"$COVERAGE"_samp_1.fq
seqtk sample -s9 ../data/6909_2.fastq.gz $READS > ../results/"$COVERAGE"_samp_2.fq

# merge and cleanup
cat ../results/"$COVERAGE"_samp_1.fq ../results/"$COVERAGE"_samp_2.fq | gzip > ../results/"$COVERAGE"_samp.fq.gz
rm ../results/"$COVERAGE"_samp_1.fq ../results/"$COVERAGE"_samp_2.fq
