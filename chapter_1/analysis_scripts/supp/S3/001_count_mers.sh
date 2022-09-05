#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="count_mers"
#SBATCH --array=5-20
#SBATCH -p koeniglab

module load jellyfish/2.2.9

## Count 5-20mers
for FILE in ../data/*.fa
do
	NAME=$(basename "$FILE" | cut -d. -f1) # basename
	
        jellyfish count -C -m "$SLURM_ARRAY_TASK_ID" -s 500M -t 4 -o "$NAME"_"$SLURM_ARRAY_TASK_ID".jf "$FILE"
        jellyfish dump -t -c "$NAME"_"$SLURM_ARRAY_TASK_ID".jf | gzip > ../results/counts/"$NAME"_"$SLURM_ARRAY_TASK_ID".txt.gz
        
        rm "$NAME"_"$SLURM_ARRAY_TASK_ID".jf
        
done
