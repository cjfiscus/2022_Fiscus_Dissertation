#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --job-name="sr"
#SBATCH -p koeniglab

conda deactivate
module unload miniconda2
module load miniconda3
source activate PyEnv3

TABLE=/rhome/dkoenig/bigdata/SCRATCH/MODEL_FIT_SEQ_BY_REMOVE/DATA/OUTPUT/FULL_MATRIX/athaliana_12mers_filtered_gc_qq_final_corrected.txt
IN=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/data/sequences/simple_repeats2.txt
OUT=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/feat_abund/simple_repeats2.txt
ABUN_TABLE=../../results/feat_abund/A_thal_kmer_abund2.txt

# export simple repeat kmers from table
python 005a_exp_simple_repeats.py "$IN" "$TABLE" > "$OUT"

# add to table
#cat "$ABUN_TABLE" <(tail -n+2 "$OUT") > ../../results/feat_abund/temp.txt
#mv "$ABUN_TABLE" ../../results/feat_abund/old3.txt
#mv ../../results/feat_abund/temp.txt "$ABUN_TABLE"
