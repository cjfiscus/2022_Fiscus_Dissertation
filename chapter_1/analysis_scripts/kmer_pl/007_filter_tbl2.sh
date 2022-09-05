#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="filter"
#SBATCH --time=1-00:15:00
#SBATCH -p batch

module load plink/1.90b3.38

KMER_TABLE="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered_gc_qq.txt"
FILTERED="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered_gc_qq_final.txt"
GC="../results/filters/gc_cor_postqq.txt"
GC_FILT="../results/filters/gc_corr_outliers.txt"
VCF="/rhome/cfisc004/shared/VCF_BIG_PANELS/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"

# aggregate gc corrs
#cat ../results/gc/parts/*.txt > "$GC" 

# identify corr outliers
#Rscript 007a_gc_cor_filter.R "$GC" "$GC_FILT"

# generate lst of remaining IDs
#head -n1 "$KMER_TABLE" | sed 's/\t/\n/g' | tail -n+2 > ../results/filters/ids.txt
#grep -v -w -f "$GC_FILT" ../results/filters/ids.txt > ../results/filters/keep_id.txt

# id of near identical lines (from genotypes)
## calc dist mtx using plink
#plink --distance square gz --vcf "$VCF" --out dist

## filter near identical lines from matrix
#Rscript 007b_id_similar_lines.R ../results/filters/

# filter table
Rscript 007c_apply_filters2.R "$KMER_TABLE" "$FILTERED"
