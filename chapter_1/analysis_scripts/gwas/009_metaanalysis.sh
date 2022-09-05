#!/bin/bash
#SBATCH --job-name=meta
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=2
#SBATCH --mem=200gb
#SBATCH -t 06:00:00
#SBATCH -p koeniglab

# environment
#conda deactivate
#module unload miniconda2
#source activate PyEnv3

LIST=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/gwas/meta/sr_lst2.txt
SCRIPT=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/scripts/gwas/009a_x2.py
SCRIPT2=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/scripts/gwas/009b_beta_mean.py
OUT=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/gwas/meta/sr_meta_cov.txt
OUT2=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/gwas/meta/sr_beta_mean_cov.txt
WD=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/gwas/assoc3/

# move into WD
cd $WD

# run scripts
python "$SCRIPT" "$LIST" > "$OUT"
python "$SCRIPT2" "$LIST" > "$OUT2"
