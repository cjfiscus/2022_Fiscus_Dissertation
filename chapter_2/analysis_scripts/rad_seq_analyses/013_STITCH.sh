#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=std/stitch%j.stdout
#SBATCH --error=std/stitch%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="impute"
#SBATCH -p koeniglab
#SBATCH --array=1-16

# software dependencies
# stitch 1.6.6
DATA=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/stitch
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stitch
CHR="SCF_$SLURM_ARRAY_TASK_ID"
K=30

# impute
STITCH.R --chr="$CHR" --bamlist="$DATA"/bamlist.txt --posfile="$DATA"/pos_"$CHR".txt --outputdir="$RESULTS" --K="$K" --nGen=100 --nCores=4
