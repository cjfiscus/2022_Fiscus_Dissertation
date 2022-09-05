#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:15:00
#SBATCH --job-name="abun"
#SBATCH -p koeniglab
#SBATCH --array=2

# python 2 environment
conda deactivate
module unload miniconda2
module load miniconda3
source activate PyEnv3

FEATURES="/rhome/cfisc004/bigdata/projects/K-mers_Arabidopsis/data/sequences/A_thal_single_copy_buscos.fa"
KMER_TABLE="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/Col-CEN_v1.2_12.txt"
RESULTS="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results"

python 002a_assign_seq_score.py \
	<(cut -f1,${SLURM_ARRAY_TASK_ID} ${KMER_TABLE} | tail -n +2) "$FEATURES" 12 \
	> "$RESULTS"/Col-CEN_v.1.2_seq_abun_buscos.txt 
