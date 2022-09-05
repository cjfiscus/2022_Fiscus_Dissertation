#!/bin/bash
#SBATCH --job-name=af
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=2
#SBATCH --mem=32gb
#SBATCH -t 06:00:00
#SBATCH -p koeniglab
#SBATCH --array=1-10

# environment
module load plink/1.90b3.38
VCF=/rhome/cfisc004/shared/VCF_BIG_PANELS/1001genomes_snp-short-indel_only_ACGTN_filtered.vcf.gz
FILE=$(ls ../../data/group_lsts/*.txt | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1)
OUT=$(basename "$FILE" | cut -d"_" -f2 | sed 's/.txt//g' | sed 's/ /_/g')

# calc af per group
plink --vcf "$VCF" --set-missing-var-ids @:# --freq --keep "$FILE" --out ../../results/af/"$OUT" 
