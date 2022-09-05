#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="stats"
#SBATCH -p koeniglab

# software dependencies
## samtools 1.10

# load required modules (slurm) 
module load samtools/1.10 

DIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/processed_bams
OUTDIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results
OUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/mapstats/cbp_mapstats.txt

## calculate flagstats
#for BAM in "$DIR"/*.bam
#do
#	BN=$(basename "$BAM")
#	samtools flagstat "$BAM" > "$OUTDIR"/mapstats/"$BN".flagstat
#done

## add to output
#for STAT in "$OUTDIR"/mapstats/*.flagstat
#do
#	python 004a_mapstat.py "$STAT" >> "$OUT"
#done

## calculate idx stats
for BAM in "$DIR"/*.bam
do
    BN=$(basename "$BAM")
    samtools idxstats "$BAM" > "$OUTDIR"/idxstats/"$BN".idxstats
done
