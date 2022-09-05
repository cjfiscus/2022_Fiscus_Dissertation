#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=std/struc%j.stdout
#SBATCH --error=std/struc%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="struc"
#SBATCH -p intel
#SBATCH --array=1-20

# software dependencies
# structure 2.3.4

module load structure/2.3.4

# define vars
DIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/structure
INPUT=input.structure
OUTPUT=output_"$SLURM_ARRAY_TASK_ID"_rep3
SEED=42069
K="$SLURM_ARRAY_TASK_ID"
N_IND=1255
N_LOC=13425

cd "$DIR"
# run structure
structure -i "$INPUT" -o "$OUTPUT" -D "$SEED" -m mainparams -K "$K" -N "$N_IND" -L "$N_LOC"
