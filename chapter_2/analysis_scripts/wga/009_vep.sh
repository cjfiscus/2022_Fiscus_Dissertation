#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="vep"
#SBATCH -p koeniglab
#SBATCH --time=4-00:00:00

conda activate vep
INPUT=/rhome/cfisc004/bigdata/projects/capsella_genomes/data/co_vep_input.txt
OUTPUT=Co_vep.txt
GFF=/rhome/cfisc004/bigdata/projects/capsella_genomes/data/Co39_temp.gff.gz
GENOME=/rhome/cfisc004/shared/CAPSELLA_GENOMES/genomes/Co39_plus_manual_correction.fasta

vep --force_overwrite -v -i "$INPUT" -o "$OUTPUT" --gff "$GFF" --fasta "$GENOME" --tab
