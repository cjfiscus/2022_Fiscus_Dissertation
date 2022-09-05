#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=std/dl%j.stdout
#SBATCH --error=std/dl%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=6-00:00:00
#SBATCH --job-name="dl"
#SBATCH -p koeniglab

cd ../../data/
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR175/000/SRR1751470/SRR1751470_2.fastq.gz
