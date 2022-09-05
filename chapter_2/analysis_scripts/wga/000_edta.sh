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
GENOME=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Cbp2-2_plus_manual_correction.fasta
EDTA.pl --genome "$GENOME" --overwrite 1 --anno 1 --sensitive 1 --anno 1 --evaluate 1 --threads 9
