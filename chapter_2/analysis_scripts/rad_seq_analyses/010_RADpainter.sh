#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=std/rad%j.stdout
#SBATCH --error=std/rad%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --job-name="rad"
#SBATCH -p koeniglab

# software dependencies
# RADpainter 0.3.3 r111
# gcc 8.3.0

module load -s base/gcc/8.3.0

INPUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/stacks/pops5/populations.haps.radpainter

# coancestry mtx
cd /rhome/cfisc004/bigdata/projects/capsella_radseq/results/RADpainter
RADpainter paint "$INPUT"
