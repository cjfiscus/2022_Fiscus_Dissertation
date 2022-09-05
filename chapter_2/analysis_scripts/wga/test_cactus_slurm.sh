#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="test_c"
#SBATCH -p koeniglab 

module unload python perl
module unload miniconda2
module load miniconda3
conda activate cactus

export TOIL_SLURM_ARGS="-p koeniglab"

cactus-prepare-toil ./cactus \
	/rhome/cfisc004/software/cactus-bin-v1.2.3/examples/evolverMammals.txt \
	--batchSystem Slurm --realTimeLogging --outHal /rhome/cfisc004/software/cactus-bin-v1.2.3/examples/evolverMammals.hal \
	--defaultDisk 20G --defaultMemory 12G --defaultCores 4 --binariesMode local
