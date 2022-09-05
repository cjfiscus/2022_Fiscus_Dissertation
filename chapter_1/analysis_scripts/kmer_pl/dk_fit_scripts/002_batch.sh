#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --array=1-840%840
#SBATCH --output=correctkmer
#SBATCH --job-name="correctkmer"
#SBATCH --partition=batch
cd /rhome/dkoenig/bigdata/SCRATCH/MODEL_FIT_SEQ_BY_REMOVE/DATA/OUTPUT/CORRECTED
thefile=`ls ../SPLIT_IT/ | head -n ${SLURM_ARRAY_TASK_ID} | tail -n1`
Rscript  ../../../SCRIPTS/R_proc.r /rhome/dkoenig/shared/MANUSCRIPT_FILES/2020_FISCUS/Data_Repository/tables/accessions/All_Samples.txt header.txt ../SPLIT_IT/${thefile} out.${thefile}
