#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=std/admix%j.stdout
#SBATCH --error=std/admix%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="admix"
#SBATCH -p koeniglab
#SBATCH --array=1-10

# software dependencies
# admixture 1.3.0

# define vars
DIR=/rhome/cfisc004/bigdata/cowpea_gwas/results/admixture
INPUT=/rhome/cfisc004/bigdata/cowpea_gwas/results/stats/pruned.bed
K="$SLURM_ARRAY_TASK_ID"

cd "$DIR"
# run admixture for 10 reps
for r in {1..10}
do
	admixture "$INPUT" "$K" -j4 --cv -s "$RANDOM" | tee log.${K}.r$r.out
	mv pruned."$K".P pruned."$K".r"$r".P
	mv pruned."$K".Q pruned."$K".r"$r".Q
done
