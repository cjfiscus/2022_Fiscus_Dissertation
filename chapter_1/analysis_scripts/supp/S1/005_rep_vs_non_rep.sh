#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="win"
#SBATCH -p koeniglab

module load bedtools/2.27.1

# environment
module unload miniconda2
module load anaconda3
source activate PyEnv2 # python 2

REFERENCE=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
INDEX=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai
ANNOTATION=../results/annotation/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.out.gff

# get repetitive seqs
bedtools merge -i "$ANNOTATION" > ../results/repetitive.bed
bedtools getfasta -fi "$REFERENCE" -bed ../results/repetitive.bed > ../results/repetitive.fa 

# get non-repetitive seqs
cut -f1,2 "$INDEX" | grep -vE "Mt|Pt" > ../results/genome.bed
bedtools complement -i "$ANNOTATION" -g ../results/genome.bed > ../results/non_repetitive.bed
bedtools getfasta -fi "$REFERENCE" -bed ../results/non_repetitive.bed > ../results/non_repetitive.fa

# calculate windowed residual 
# variables
SEQS=../results/repetitive.fa
TABLE=../results/residuals.txt

# calculate median residual for repetitive seq
python 005a_assign_seq_score.py <(cut -f1,2 ${TABLE} | tail -n+2) "$SEQS" 12 > ../results/residuals_repetitive.txt 

#####

# variables
SEQS=../results/non_repetitive.fa
TABLE=../results/residuals.txt

# calculate median residual for non-repetitive seq
python 005a_assign_seq_score.py <(cut -f1,2 ${TABLE} | tail -n+2) "$SEQS" 12 > ../results/residuals_non_repetitive.txt 
