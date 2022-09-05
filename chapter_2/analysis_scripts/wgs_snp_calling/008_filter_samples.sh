#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=5-00:00:00
#SBATCH --job-name="filter"
#SBATCH -p koeniglab

# software dependencies
## bcftools/1.10
## plink 1.9
module load bcftools/1.10

#### PIPELINE #####
# filter the following individuals
# Cr: 104.12
SPP=Cr
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
#bcftools view -s ^104.12 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
#tabix "$RESULTS"/"$SPP"_filtered.vcf.gz

# Cr_old: 104.12
SPP=Cr_old
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
#bcftools view -s ^104.12 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
#tabix "$RESULTS"/"$SPP"_filtered.vcf.gz

# CbpCr: 104.12
SPP=CbpCr
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
#bcftools view -s ^104.12 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
#tabix "$RESULTS"/"$SPP"_filtered.vcf.gz

# Cbp: BEL5, KYRG-3-14
SPP=Cbp
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
#bcftools view -s ^BEL5,KYRG-3-14 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
#tabix "$RESULTS"/"$SPP"_filtered.vcf.gz

# Co: 1981-10
SPP=Co
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
#bcftools view -s ^1981-10 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
#tabix "$RESULTS"/"$SPP"_filtered.vcf.gz

# CbpCo: 1981-10
SPP=CbpCo
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
#bcftools view -s ^1981-10 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
#tabix "$RESULTS"/"$SPP"_filtered.vcf.gz

# CbpCo_Cbp: 1981-10, BEL5, KYRG-3-14
SPP=CbpCo_Cbp
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
cp "$RESULTS"/"$SPP"_filter3.vcf.gz "$RESULTS"/"$SPP"_filtered.vcf.gz
#bcftools view -s ^1981-10,BEL5,KYRG-3-14 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
tabix "$RESULTS"/"$SPP"_filtered.vcf.gz

# CbpCr_Cbp: 104.12, BEL5, KYRG-3-14
SPP=CbpCr_Cbp
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
cp "$RESULTS"/"$SPP"_filter3.vcf.gz "$RESULTS"/"$SPP"_filtered.vcf.gz 
#bcftools view -s ^104.12,BEL5,KYRG-3-14 "$RESULTS"/"$SPP"_filter3.vcf.gz | bgzip > "$RESULTS"/"$SPP"_filtered.vcf.gz
tabix "$RESULTS"/"$SPP"_filtered.vcf.gz
