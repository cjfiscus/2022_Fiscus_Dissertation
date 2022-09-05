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
DIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/admixture
INPUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops4/plink.bed
K="$SLURM_ARRAY_TASK_ID"

cd "$DIR"
# run admixture for 10 reps
for r in {1..10}
do
	admixture "$INPUT" "$K" -j4 --cv -s "$RANDOM" | tee log.${K}.r$r.out
	mv plink."$K".P plink."$K".r"$r".P
	mv plink."$K".Q plink."$K".r"$r".Q
done
