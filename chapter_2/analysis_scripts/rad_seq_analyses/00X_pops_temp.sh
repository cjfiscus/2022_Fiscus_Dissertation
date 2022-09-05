#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=std/popss%j.stdout
#SBATCH --error=std/popss%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --job-name="stacks"
#SBATCH -p koeniglab

# software dependencies
# stacks 2.60

module load stacks/2.60

# define vars
THREADS=4
DIR=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks
OUTDIR1=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops3
OUTDIR2=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops4
OUTDIR3=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops5
POPMAP=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/popmap_filter2.txt

# pops call with filters
## snp based genotypes
#populations -t "$THREADS" --popmap "$POPMAP" -P "$DIR" -O "$OUTDIR1" \
#	-r 0.60 --min-maf 0.01 --max-obs-het 0.05 \
#	--vcf --ordered-export --genepop --plink --verbose

## pruned snp based genotypes (for structure)
populations -t "$THREADS" --popmap "$POPMAP" -P "$DIR" -O "$OUTDIR2" \
        -r 0.60 --min-maf 0.01 --max-obs-het 0.05 \
        --write-single-snp --ordered-export --plink --genepop --verbose

### prep structure file
#tail -n+2 "$OUTDIR2"/populations.structure > ../results/structure/input.structure

## haplotype based genotypes (for RADpainter)
#populations -t "$THREADS" --popmap "$POPMAP" -P "$DIR" -O "$OUTDIR3" \
#	-r 0.60 --min-maf 0.01 --max-obs-het 0.05 \
#	--filter-haplotype-wise --radpainter --verbose
