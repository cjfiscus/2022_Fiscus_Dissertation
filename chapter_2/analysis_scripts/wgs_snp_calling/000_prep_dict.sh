#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=16G
#SBATCH --output=std/dict%j.stdout
#SBATCH --error=std/dict%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=4-00:00:00
#SBATCH --job-name="dict"
#SBATCH -p koeniglab

# software dependencies
## GATK 4.1.8.0
## bwa-mem2 v. 2.0
## samtools/1.10

module load samtools/1.10

# SET VARIABLES
REFERENCE=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Cbp2-2_plus_manual_correction.fasta
OUT=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Cbp2-2_plus_manual_correction.dict

# create .dict
#gatk --java-options "-Xmx8G" CreateSequenceDictionary \
#            -R "$REFERENCE" \
#            -O "$OUT"

# index with bwa-mem2
#bwa-mem2 index "$REFERENCE"

# index with samtools
samtools faidx "$REFERENCE"
