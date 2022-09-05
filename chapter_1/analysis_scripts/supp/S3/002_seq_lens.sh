#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="sum"
#SBATCH -p koeniglab

module load samtools/1.9
module unload miniconda2
module load anaconda3
source activate PyEnv3

for i in ../data/*.fa
do
	samtools faidx $i
	python 002a_faidx_sum.py "$i".fai >> ../results/seq_lens.txt

done
