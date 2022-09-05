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
INPUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/gwas/plink.vcf
OUTPUT=Cbp_vep.txt
GFF=/rhome/cfisc004/bigdata/projects/capsella_genomes/data/Cbp2-2_temp.gff.gz
GENOME=/rhome/cfisc004/shared/CAPSELLA_GENOMES/genomes/Cbp2-2_plus_manual_correction.fasta

vep --force_overwrite -v -i "$INPUT" -o "$OUTPUT" --gff "$GFF" --fasta "$GENOME" --tab
