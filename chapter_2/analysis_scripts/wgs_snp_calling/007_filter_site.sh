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
module load plink/1.90b3.38
module load vcftools/0.1.15

## define vars
SPP=CbpCr_Cbp
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
THREADS=4

#### PIPELINE #####
# initial filters
## GATK best practices hard filters
bcftools filter -e'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0' \
	"$RESULTS"/"$SPP"_raw.vcf \
	> "$RESULTS"/temp.vcf
bcftools view -f.,PASS "$RESULTS"/temp.vcf \
	| bgzip > "$RESULTS"/"$SPP"_filter1.vcf.gz
rm "$RESULTS"/temp.vcf

## require 3 reads to call and keep only biallelic sites
bcftools filter -e'FMT/DP<3' -S . "$RESULTS"/"$SPP"_filter1.vcf.gz \
	| bcftools view -i 'F_MISSING<1' -m2 -M2 > "$RESULTS"/temp.vcf
##########
# filter sites with > 5% het & > 5% missing data
## calculate proportion het per site with plink
if [[ "$SPP" == "CbpCr" ]] || [[ "$SPP" == "Cr" ]] || [[ "$SPP" == "Cr_old" ]] || [[ "$SPP" == "CbpCr_Cbp" ]] || [[ "$SPP" == "CbpCo_Cbp" ]]
then
# exclude grandiflora from het calcs
	plink --vcf "$RESULTS"/temp.vcf --remove ../../data/Cg.txt \
		--freqx --allow-extra-chr --out "$RESULTS"/"$SPP"
else
	plink --vcf "$RESULTS"/temp.vcf --freqx --allow-extra-chr \
		--out "$RESULTS"/"$SPP"
fi

## list of site ids for following Rscript
bcftools query -f'%CHROM %POS\n' "$RESULTS"/temp.vcf > "$RESULTS"/sites.txt

## plot AFs and produce list of sites to filter
Rscript 007a_plot_frqx.R "$RESULTS"/"$SPP".frqx "$RESULTS"/sites.txt 0.05 "$RESULTS"/"$SPP"

## set of sites to remove with het > 5%
cut -f1,2 "$RESULTS"/"$SPP"_hetmin.txt | tail -n+2 > "$RESULTS"/remove.txt

## apply filter
vcftools --vcf "$RESULTS"/temp.vcf \
	--exclude-positions "$RESULTS"/remove.txt \
	--recode --recode-INFO-all --stdout \
	| bcftools view -i 'F_MISSING<0.05' \
        | bgzip > "$RESULTS"/"$SPP"_filter2.vcf.gz

##########
# top cut depth at value of > Q3 + 1.5IQR
## calculate depth per site and plot
bcftools query -f '%CHROM %POS %DP\n' "$RESULTS"/"$SPP"_filter2.vcf.gz \
	> "$RESULTS"/depth.txt
Rscript 007b_plot_dp.R "$RESULTS"/depth.txt "$RESULTS"/"$SPP"

## filter on depth
MAXDP=$(head -n1 "$RESULTS"/"$SPP"_depth_cutval.txt)
bcftools view -i "INFO/DP < $MAXDP" "$RESULTS"/"$SPP"_filter2.vcf.gz \
	| bgzip > "$RESULTS"/"$SPP"_filter3.vcf.gz

##########
# write filtering report
touch "$RESULTS"/log.txt
echo "$SPP" >> "$RESULTS"/log.txt
echo "GATK best practices filter" >> "$RESULTS"/log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$RESULTS"/"$SPP"_filter1.vcf.gz | wc -l) \
	>> "$RESULTS"/log.txt

echo "biallelic sites, het, and missing data filter" >> "$RESULTS"/log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$RESULTS"/"$SPP"_filter2.vcf.gz | wc -l) \
	>> "$RESULTS"/log.txt

echo "depth filter " >> "$RESULTS"/log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$RESULTS"/"$SPP"_filter3.vcf.gz | wc -l) \
	>> "$RESULTS"/log.txt

##########
# calculate stats per ind
vcftools --gzvcf "$RESULTS"/"$SPP"_filter3.vcf.gz --het --stdout \
	> "$RESULTS"/"$SPP"_het.txt
Rscript 007c_plot_het.R "$RESULTS"/"$SPP"_het.txt "$RESULTS"/"$SPP"
