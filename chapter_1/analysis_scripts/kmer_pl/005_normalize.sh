#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=400gb
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --time=4-00:00:00
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="normrle"
#SBATCH -p koeniglab

# count tables 
FILTERED="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered.txt"
GCNORM="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered_gc.txt"
GCTMMNORM="/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers_filtered_gc_rle.txt"
REF_PROP_TABLE="../results/filters/ref_gc_prop.txt"
PROP_TABLE="../results/filters/raw_gc_prop.txt"

# gc correct counts
#echo "GC correct"
#Rscript 005a_gc_correct.R "$FILTERED" "$GCNORM" "$REF_PROP_TABLE" "$PROP_TABLE"  

# qq normalize counts
echo "qq norm"
Rscript 005b_rle.R "$GCNORM" "$GCTMMNORM"
