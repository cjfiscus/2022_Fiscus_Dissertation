#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
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
make_panTElib.pl -liblist liblst.txt -threads 3
