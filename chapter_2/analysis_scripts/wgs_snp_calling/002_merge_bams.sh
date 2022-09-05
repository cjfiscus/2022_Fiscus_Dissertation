#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G
#SBATCH --output=std/alignments%j.stdout
#SBATCH --error=std/alignments%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="merge"
#SBATCH -p koeniglab

# software dependencies
module load samtools/1.10

# merge bams for samples with multiple seq runs
## cbp
#samtools merge ../../results/alignments/Cbp/wt-msu.bam \
#	../../results/alignments/Cbp/SRR5412134.bam \
#	../../results/alignments/Cbp/SRR5412135.bam \
#	../../results/alignments/Cbp/SRR5412136.bam \
#	../../results/alignments/Cbp/SRR5412137.bam \
#	../../results/alignments/Cbp/SRR5803007.bam \
#	../../results/alignments/Cbp/SRR5803008.bam \
#	../../results/alignments/Cbp/SRR5803009.bam 

## CbpCr
#samtools merge ../../results/alignments/CbpCr/Monte_Gargano.bam \
#	../../results/alignments/CbpCr/SRR207609.bam \
#	../../results/alignments/CbpCr/SRR207610.bam \
#	../../results/alignments/CbpCr/SRR207611.bam \
#	../../results/alignments/CbpCr/SRR207612.bam \
#	../../results/alignments/CbpCr/SRR207619.bam \
#	../../results/alignments/CbpCr/SRR207620.bam \
#	../../results/alignments/CbpCr/SRR207621.bam \
#	../../results/alignments/CbpCr/SRR207622.bam \
#	../../results/alignments/CbpCr/SRR207623.bam \
#	../../results/alignments/CbpCr/SRR207624.bam \
#	../../results/alignments/CbpCr/SRR207625.bam \
#	../../results/alignments/CbpCr/SRR207626.bam

## Cr
#samtools merge ../../results/alignments/Cr/Monte_Gargano.bam \
#	../../results/alignments/Cr/SRR207609.bam \
#	../../results/alignments/Cr/SRR207610.bam \
#	../../results/alignments/Cr/SRR207611.bam \
#	../../results/alignments/Cr/SRR207612.bam \
#	../../results/alignments/Cr/SRR207619.bam \
#	../../results/alignments/Cr/SRR207620.bam \
#	../../results/alignments/Cr/SRR207621.bam \
#	../../results/alignments/Cr/SRR207622.bam \
#	../../results/alignments/Cr/SRR207623.bam \
#	../../results/alignments/Cr/SRR207624.bam \
#	../../results/alignments/Cr/SRR207625.bam \
#	../../results/alignments/Cr/SRR207626.bam

## Cr_old
#samtools merge ../../results/alignments/Cr_old/Monte_Gargano.bam \
#	../../results/alignments/Cr_old/SRR207609.bam \
#	../../results/alignments/Cr_old/SRR207610.bam \
#	../../results/alignments/Cr_old/SRR207611.bam \
#	../../results/alignments/Cr_old/SRR207612.bam \
#	../../results/alignments/Cr_old/SRR207619.bam \
#	../../results/alignments/Cr_old/SRR207620.bam \
#	../../results/alignments/Cr_old/SRR207621.bam \
#	../../results/alignments/Cr_old/SRR207622.bam \
#	../../results/alignments/Cr_old/SRR207623.bam \
#	../../results/alignments/Cr_old/SRR207624.bam \
#	../../results/alignments/Cr_old/SRR207625.bam \
#	../../results/alignments/Cr_old/SRR207626.bam

samtools merge ../../results/alignments/CbpCo/Monte_Gargano.bam \
        ../../results/alignments/CbpCo/SRR207609.bam \
        ../../results/alignments/CbpCo/SRR207610.bam \
        ../../results/alignments/CbpCo/SRR207611.bam \
        ../../results/alignments/CbpCo/SRR207612.bam \
        ../../results/alignments/CbpCo/SRR207619.bam \
        ../../results/alignments/CbpCo/SRR207620.bam \
        ../../results/alignments/CbpCo/SRR207621.bam \
        ../../results/alignments/CbpCo/SRR207622.bam \
        ../../results/alignments/CbpCo/SRR207623.bam \
        ../../results/alignments/CbpCo/SRR207624.bam \
        ../../results/alignments/CbpCo/SRR207625.bam \
        ../../results/alignments/CbpCo/SRR207626.bam

