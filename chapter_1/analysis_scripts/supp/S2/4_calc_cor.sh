#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=200GB
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="cor"
#SBATCH -p koeniglab

Rscript 4a_calc_cor.R ../results/tables/12mers.txt

#for i in ../results/tables/*mers.txt 
#do
#	Rscript 4a_calc_cor.R $i
#done
