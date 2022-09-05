#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200gb
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="combine"
#SBATCH -p koeniglab

## sort file by K-mer 
sort <(zcat ../results/TAIR10.txt.gz) > ../results/TAIR10.tmp
sort <(zcat ../results/6909.txt.gz) > ../results/6909.tmp

## merge first two files to create master
join -1 1 -2 1 -a 1 -a 2 -e 0 -o 0,1.2,2.2 -t $'\t' ../results/TAIR10.tmp ../results/6909.tmp \
	> ../results/"$SLURM_JOB_ID".tmp

## add header
HEADER=$(echo "mer\tTAIR10\t6909")
sed '1s/^/'"$HEADER"'\n/' ../results/"$SLURM_JOB_ID".tmp > ../results/"$SLURM_JOB_ID"_2.tmp
mv ../results/"$SLURM_JOB_ID"_2.tmp ../results/"$SLURM_ARRAY_TASK_ID"mers.txt

# cleanup
rm ../results/*.tmp
