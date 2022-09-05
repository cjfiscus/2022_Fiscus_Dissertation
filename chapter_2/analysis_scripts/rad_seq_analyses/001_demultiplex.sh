#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="demux"
#SBATCH -p koeniglab
#SBATCH --array=2-7

# software dependencies
# stacks 2.60

module load stacks/2.60

## ARRAY 2-7

#VARS
LST=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/manifest.txt
FILE_IN=$(head -n "$SLURM_ARRAY_TASK_ID" "$LST" | tail -n 1 | cut -f3)
DIR_OUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/demulti/$(echo "$FILE_IN" | awk -F/ '{print $(NF-1)}')
BARCODE=$(head -n "$SLURM_ARRAY_TASK_ID" "$LST" | tail -n 1 | cut -f2)
ENZ=kpnI

# make output folder
mkdir "$DIR_OUT"

# demultiplex 
# -r rescue barcodes, -c rm reads with uncalled base, -q discard reads with low qual
process_radtags -f "$FILE_IN" -o "$DIR_OUT" -b "$BARCODE" -e "$ENZ" -r -c -q
