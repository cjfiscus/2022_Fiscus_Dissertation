#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="stat"
#SBATCH -p koeniglab
#SBATCH --array=5-20

module unload miniconda2
module load anaconda3
source activate PyEnv3

# calculate stats
python 003a_stats.py ../results/counts/*_"$SLURM_ARRAY_TASK_ID".txt.gz > ../results/stats/"$SLURM_ARRAY_TASK_ID".txt
