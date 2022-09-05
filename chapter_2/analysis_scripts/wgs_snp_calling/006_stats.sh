#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=5-00:00:00
#SBATCH --job-name="plotter"
#SBATCH -p koeniglab

# software dependencies
## bcftools/1.10
module load bcftools/1.10

## vars
SPP=CbpCr_Cbp
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
THREADS=4

#### PIPELINE #####
# calculate stats
bcftools query -f '%CHROM %POS %REF %ALT %QD %FS %SOR %MQ %MQRankSum %ReadPosRankSum\n' \
	"$RESULTS"/"$SPP"_raw.vcf | gzip > "$RESULTS"/"$SPP"_raw_stats.txt.gz

# plot 
Rscript 006a_plot_stats.R "$RESULTS"/"$SPP"_raw_stats.txt.gz 0.10 "$RESULTS"/"$SPP"
