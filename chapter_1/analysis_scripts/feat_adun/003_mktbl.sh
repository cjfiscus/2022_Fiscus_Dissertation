#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=01:00:00
#SBATCH --job-name="short"
#SBATCH -p short

conda deactivate
module unload miniconda2
module load miniconda3
source activate PyEnv3

KMER_TABLE="/rhome/dkoenig/bigdata/SCRATCH/MODEL_FIT_SEQ_BY_REMOVE/DATA/OUTPUT/FULL_MATRIX/athaliana_12mers_filtered_gc_qq_final_corrected.txt"
HEADER="../../results/feat_abund/header.txt"
TAILER="../../results/feat_abund/tailer.txt"
OUT="../../results/feat_abund/A_thal_kmer_abund.txt"

# make header 
head -n 1 "$KMER_TABLE" | cut -f2- | sed '1s/^/Feature\t/' > "$HEADER" 

# combine files in column-wise manner (in naming order) 
python 003a_merge.py > "$TAILER" 

# add feature names
paste <(cut -f1 ../../results/feat_abund/parts/part_2.txt) "$TAILER" > "$TAILER".tmp 

# join header and tailer 
cat "$HEADER" "$TAILER".tmp  > "$OUT"
rm "$TAILER".tmp
