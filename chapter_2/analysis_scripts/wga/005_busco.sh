#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="busco"
#SBATCH --time=3-00:00:00
#SBATCH -p koeniglab

# set environment
conda deactivate
module load busco/5.3.2

THREADS=4
GENOME=/rhome/cfisc004/bigdata/projects/capsella_genomes/scripts/wga/Capsella_orientalis.fa
LINEAGE=/rhome/cfisc004/bigdata/projects/capsella_genomes/scripts/wga/busco_downloads/lineages/brassicales_odb10
OUT=Co_Agren
OUTDIR=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/busco

busco -c 4 -m genome -i "$GENOME" -l "$LINEAGE" -o "$OUT" --out_path "$OUTDIR" 
