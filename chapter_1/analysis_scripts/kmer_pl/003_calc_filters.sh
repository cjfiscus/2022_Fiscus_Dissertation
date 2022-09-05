#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100G
#SBATCH --output=./std/pl%j.stdout
#SBATCH --error=./std/pl%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="filter"
#SBATCH --time=2-00:00:00
#SBATCH -p koeniglab

conda deactivate

module unload miniconda2
module load anaconda3
source activate PyEnv3

KMERTABLE=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers.txt

## mapping % filter
# combine mapping stats into one file 
for F in ../results/mappingstats/genome/*.txt
do
        python ./003a_mapstat.py $F >> ../results/mappingstats/mapstats_genome.txt
done
## remove NAs
sed 's/N\/A/0/g' ../results/mappingstats/mapstats_genome.txt > temp.txt
mv temp.txt ../results/mappingstats/mapstats_genome.txt

# plot and write out filter set
Rscript ./003b_mapping_filter.R ../results/mappingstats/mapstats_genome.txt ../results/filters/low_mapping_ref.txt

## coverage filter
Rscript ./003c_coverage_filter.R "$KMERTABLE" ../results/filters/coverage.txt ../results/filters/low_est_coverage.txt

## extreme gc filter
# calculate gc
python ./003d_calc_gc_lib.py "$KMERTABLE" ../results/filters/gc_content_lib.txt 

# filter on gc
Rscript 003e_gc_filter.R ../results/filters/gc_content_lib.txt ../results/filters/gc_extremes.txt

# cat initial filter set 
cd ../results/filters
cat low_est_coverage.txt gc_extremes.txt low_mapping_ref.txt | grep -v "lyrata" | grep -v "rubella" | sort -u > filtered1.txt
