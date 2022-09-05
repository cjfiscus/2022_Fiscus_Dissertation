#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=100G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="filter_vcf"
#SBATCH -p koeniglab 

module load bcftools/1.8

VCF_IN=/rhome/cfisc004/shared/VCF_BIG_PANELS/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
VCF_OUT=/rhome/cfisc004/shared/VCF_BIG_PANELS/1001genomes_snp-short-indel_only_ACGTN_filtered.vcf.gz
ANNOT_TBL=../../results/feat_abund/A_thal_kmer_abund.txt
PASS_FILT=../../results/gwas/id_lst.txt

# extract list of samples to keep 
## parse ids from seq abun table
head -n 1 "$ANNOT_TBL" | tr '\t' '\n' | tail -n+2 | sort -n > ../../results/gwas/temp1.txt

## parse ids from vcf
bcftools query -l "$VCF_IN" > ../../results/gwas/temp2.txt

## find accessions in common
join <(sort ../../results/gwas/temp1.txt) <(sort ../../results/gwas/temp2.txt) | sort -n > "$PASS_FILT"

# Filter VCF as follows: 
## subset samples that passed filtering (--samples-file "$PASS_FILT")
## no organellar calls (--regions 1,2,3,4,5)
## minimum quality of 30 (-i 'MIN(QUAL)>30')
## only biallelic variants [min 2, max 2] (-m2 -M2)

bcftools view --samples-file "$PASS_FILT" --regions 1,2,3,4,5 -i 'MIN(QUAL)>30'\
	-m2 -M2 "$VCF_IN" |\
	bgzip > "$VCF_OUT"

## index vcf 
tabix "$VCF_OUT" 
