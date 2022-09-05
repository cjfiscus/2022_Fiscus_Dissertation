#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=64G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="verify"
#SBATCH --time=3-00:00:00
#SBATCH -p koeniglab

# set environment
module unload python perl
module unload miniconda2
module load miniconda3
conda activate cactus

## variables
HAL=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/alignments/msa/aln5.hal

# validate
halValidate $HAL

# calc high level stats
halStats $HAL > $HAL.stats

# summarize mutations
halSummarizeMutations $HAL > $HAL.mut

# convert to maf
hal2maf $HAL $HAL.maf
