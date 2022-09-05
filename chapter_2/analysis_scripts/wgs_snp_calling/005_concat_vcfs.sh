#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G
#SBATCH --output=std/cat%j.stdout
#SBATCH --error=std/cat%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=5-00:00:00
#SBATCH --job-name="concat"
#SBATCH -p koeniglab

# software dependencies
## bcftools/1.10
module load bcftools/1.10

## vars
SPP=CbpCo_Cbp
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
THREADS=4

#### PIPELINE #####

# concatonate VCFs
bcftools concat --threads 3 -f file_lst_"$SPP".txt -o "$RESULTS"/"$SPP"_raw.vcf 
