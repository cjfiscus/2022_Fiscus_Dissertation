#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="count_mers"
#SBATCH -p koeniglab
#SBATCH --array=12

module load jellyfish/2.2.9

REFERENCE=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

cd ../results/

# count in reference genome
jellyfish count -C -m "$SLURM_ARRAY_TASK_ID" -s 500M -t 4 -o TAIR10_"$SLURM_ARRAY_TASK_ID".jf "$REFERENCE"
jellyfish dump -t -c TAIR10_"$SLURM_ARRAY_TASK_ID".jf | gzip > ./counts/"$SLURM_ARRAY_TASK_ID"/TAIR10.txt.gz
rm TAIR10_"$SLURM_ARRAY_TASK_ID".jf

# count in simulated or sampled reads
for FILE in *.fq.gz 
do
	NAME=$(basename "$FILE" | cut -d. -f1)		
	jellyfish count -C -m "$SLURM_ARRAY_TASK_ID" -s 500M -t 4 -o "$NAME"_"$SLURM_ARRAY_TASK_ID".jf <(zcat "$FILE") 
        jellyfish dump -t -c "$NAME"_"$SLURM_ARRAY_TASK_ID".jf | gzip > ./counts/"$SLURM_ARRAY_TASK_ID"/"$NAME".txt.gz
       
	# cleanup 
        rm "$NAME"_"$SLURM_ARRAY_TASK_ID".jf
done
