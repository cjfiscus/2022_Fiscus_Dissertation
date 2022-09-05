#!/bin/zsh -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --job-name="ext"
#SBATCH --partition=koeniglab

module load samtools/1.9

# python 3
conda deactivate
module unload miniconda2
module load anaconda3
source activate PyEnv3

GENOME=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
INDEX="$GENOME".fai
WINDOW=100000
OUT=../../data/sequences/TAIR10_"$WINDOW".txt
OUT_FA=../../data/sequences/TAIR10_"$WINDOW".fa

# generate intervals to extract from genome
python 001a_gen_intervals.py "$INDEX" "$WINDOW" | grep -v "Mt\|Pt" > "$OUT"

# extract sequences
xargs samtools faidx "$GENOME" < "$OUT" > "$OUT_FA"
