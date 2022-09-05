#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=400G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="pca"
#SBATCH --time=1-00:15:00
#SBATCH -p koeniglab

TABLE="/rhome/dkoenig/bigdata/SCRATCH/MODEL_FIT_SEQ_BY_REMOVE/DATA/OUTPUT/FULL_MATRIX/athaliana_12mers_filtered_gc_qq_final_corrected.txt"

# do pca
Rscript 008a_pca.R "$TABLE" qq_corrected
