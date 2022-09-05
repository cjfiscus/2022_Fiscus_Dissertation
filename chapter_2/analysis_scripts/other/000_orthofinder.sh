#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="ortho"
#SBATCH -p koeniglab

# import module
module load orthofinder/2.5.4
IN_DIR=/rhome/cfisc004/bigdata/projects/capsella_genomes/data/proteomes
OUT_DIR=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/orthofinder

orthofinder -t 16 -a 16 -f "$IN_DIR" -o "$OUT_DIR"
