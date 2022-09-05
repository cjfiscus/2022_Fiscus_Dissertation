#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=std/popss%j.stdout
#SBATCH --error=std/popss%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --job-name="stacks"
#SBATCH -p koeniglab

# software dependencies
# stacks 2.60

module load stacks/2.60

# define vars
THREADS=4
DIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks
OUTDIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops2
POPMAP=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/popmap_filter1.txt

# pops call with filters
populations -t "$THREADS" --popmap "$POPMAP" -P "$DIR" -O "$OUTDIR" -r 0.60 --min-maf 0.01 --max-obs-het 0.05 --vcf --verbose
