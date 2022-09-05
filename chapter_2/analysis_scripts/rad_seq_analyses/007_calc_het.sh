#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=64G
#SBATCH --output=std/vcftools%j.stdout
#SBATCH --error=std/vcftools%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --job-name="vcftools"
#SBATCH -p koeniglab

# software dependencies
# vcftools 0.1.15

module load vcftools/0.1.15

# define vars
VCF=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops1/populations.snps.vcf
OUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops1/het_by_sample.txt

# calculate het calls with vcftools
vcftools --vcf "$VCF" --het --stdout > "$OUT"
