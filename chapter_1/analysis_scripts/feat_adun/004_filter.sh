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

TABLE="../../results/feat_abund/A_thal_kmer_abund.txt"
RESULT_DIR="../../results/feat_abund"

## filter non-representative non-coding genes
Rscript 004a_filter.R "$TABLE" "$RESULT_DIR"/A_thal_kmer_abund_filtered.txt

# AT1G06740.1 is part of protein-coding genes and TE genes so extra value needs to be removed
# AT3G05850.1 is part of protein-coding genes and TE genes so extra value needs to be removed
# remove mitochondrial genes
awk '/^AT1G06740\.1/ && !f{f=1; next} 1' "$RESULT_DIR"/A_thal_kmer_abund_filtered.txt |
awk '/^AT3G05850\.1/ && !f{f=1; next} 1' | grep -v "^ATMG" > "$RESULT_DIR"/temp.txt

# move unfiltered table to old, filtered table to original name
mv "$TABLE" "$RESULT_DIR"/old.txt
mv "$RESULT_DIR"/temp.txt "$TABLE"
rm "$RESULT_DIR"/A_thal_kmer_abund_filtered.txt
