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

# environment
module load plink/1.90b3.38
VCF=/rhome/cfisc004/shared/VCF_BIG_PANELS/1001genomes_snp-short-indel_only_ACGTN_filtered.vcf.gz

# calc all af
plink --vcf "$VCF" --set-missing-var-ids @:# --freq --out ../../results/af/all
