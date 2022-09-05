#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=pl%j.stdout
#SBATCH --error=pl%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="count"
#SBATCH -p batch

# software dependencies
## axel 2.16.1 (https://github.com/axel-download-accelerator/axel/releases) 
## jellyfish 2.2.9; 

# load required modules (slurm) 
module load jellyfish/2.2.9

# SET VARIABLES
WORKINGDIR=../
RESULTS_DIR=/rhome/cfisc004/bigdata/projects/SX/results/
TEMP_DIR=../data/
THREADS=4

#### PIPELINE #####

# make temp directory and go there
TEMP_DIR="$TEMP_DIR""$NAME"
mkdir -pv "$TEMP_DIR"
cd "$TEMP_DIR"

## download assemblies
axel -n "$THREADS" http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/om/OMOK01.fasta.gz -o Sequel_assem.fastq.gz

axel -n "$THREADS" ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/om/OMOL01.fasta.gz -o MinION_assem.fastq.gz

# count k-mers
jellyfish count -C -m 12 -s 500M -t 4 -o "$RESULTS_DIR"/Sequel_assem.jf <(zcat Sequel_assem.fastq.gz)
jellyfish dump -tc "$RESULTS_DIR"/Sequel_assem.jf | gzip > "$RESULTS_DIR"/Sequel_assem.txt.gz

jellyfish count -C -m 12 -s 500M -t 4 -o "$RESULTS_DIR"/MinION_assem.jf <(zcat MinION_assem.fastq.gz)
jellyfish dump -tc "$RESULTS_DIR"/MinION_assem.jf | gzip > "$RESULTS_DIR"/MinION_assem.txt.gz

rm "$RESULTS_DIR"/*.jf
