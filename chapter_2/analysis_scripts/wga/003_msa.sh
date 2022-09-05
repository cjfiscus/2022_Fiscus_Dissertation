#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="msa1"
#SBATCH --time=14-00:00:00
#SBATCH -p koeniglab

# set environment
module unload python perl
module unload miniconda2
module load miniconda3
conda activate cactus

## variables
VAR=1
JOBSTORE=./cactus_"$VAR"
SEQFILE=/rhome/cfisc004/bigdata/projects/capsella_genomes/data/msa"$VAR".txt
OUT=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/alignments/msa/aln"$VAR".hal
TMPDIR=/rhome/cfisc004/bigdata/projects/capsella_genomes/scripts/wga/work/tmp
export TOIL_SLURM_ARGS="-p koeniglab --time=5-00:00:00"
export TMPDIR
# --workDir /scratch/$USER/$SLURM_JOB_ID

# align
cactus-prepare-toil --restart "$JOBSTORE" "$SEQFILE" --outHal "$OUT" \
	--batchSystem Slurm --realTimeLogging --disableCaching \
	--defaultDisk 150G --workDir /rhome/cfisc004/bigdata/projects/capsella_genomes/scripts/wga/work \
	--defaultMemory 32G --defaultCores 4 --binariesMode local
