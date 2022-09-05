#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="edta"
#SBATCH -p koeniglab
#SBATCH --time=4-00:00:00

conda deactivate
module unload python perl
module unload miniconda2
module load miniconda3
conda activate EDTA

## EDTA 2.0
THREADS=9
LIB=/rhome/cfisc004/bigdata/projects/capsella_genomes/scripts/wga/liblst.txt.panTE.fa
GENOME=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Cbp2-2Cr_plus_manual_correction.fasta

# mask with RepeatMasker
RepeatMasker -q -div 40 -lib "$LIB" -gff "$GENOME"

# final annot with EDTA
EDTA.pl --t "$THREADS" --genome "$GENOME" --step final --anno 1 --rmout "$GENOME".out --curatedlib "$LIB"
