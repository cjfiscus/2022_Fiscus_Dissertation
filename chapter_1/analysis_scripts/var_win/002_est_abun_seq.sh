#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=16G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00
#SBATCH --job-name="abun"
#SBATCH -p koeniglab
#SBATCH --array=2-1048

# python 2 environment
conda deactivate
module unload miniconda2
module load miniconda3
source activate PyEnv3

FEATURES="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/data/sequences/TAIR10_100000.fa"
KMER_TABLE="/rhome/dkoenig/bigdata/SCRATCH/MODEL_FIT_SEQ_BY_REMOVE/DATA/OUTPUT/FULL_MATRIX/athaliana_12mers_filtered_gc_qq_final_corrected.txt"
RESULTS="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/var_win/100/"
	
python 002c_assign_sequence_score_pipe.py \
	<(cut -f1,${SLURM_ARRAY_TASK_ID} ${KMER_TABLE} | tail -n +2) "$FEATURES" 12 \
	> "$RESULTS"/part_${SLURM_ARRAY_TASK_ID}.txt
