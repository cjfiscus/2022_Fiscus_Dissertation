#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="annot"
#SBATCH -p koeniglab

# python environment for py3
module unload miniconda2
module load anaconda3
source activate PyEnv3

RAWCOUNTS=../results/12mers.txt
ANNOTCOUNTS=../results/12mers_annot.txt

# calculate GC content for libs
python 006a_calc_gc_lib.py \
	"$RAWCOUNTS" \
	../results/gc.txt

# annot table with GC content
paste <(python 006b_annot_mers_gc.py <(cut -f1 "$RAWCOUNTS") | cut -f1) "$RAWCOUNTS" > "$ANNOTCOUNTS"

# calculate correlation between GC and count
Rscript 006c_gc_cor.R 
