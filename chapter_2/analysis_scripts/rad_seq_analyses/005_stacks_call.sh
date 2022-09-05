#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100G
#SBATCH --output=std/stacks%j.stdout
#SBATCH --error=std/stacks%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=10-00:00:00
#SBATCH --job-name="stacks"
#SBATCH -p koeniglab

# software dependencies
# stacks 2.60

module load stacks/2.60

THREADS=7
DIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/processed_bams
POPMAP=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/popmap.txt
OUTDIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks
#####
ref_map.pl -T "$THREADS" --samples "$DIR" --popmap "$POPMAP" -o "$OUTDIR"
