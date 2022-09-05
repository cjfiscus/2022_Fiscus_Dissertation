#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=./std/Cbp%j.stdout
#SBATCH --error=./std/Cbp%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=5-00:00:00
#SBATCH --job-name="vcf"
#SBATCH -p koeniglab
#SBATCH --array=17

# software dependencies
# gatk 4.1.8.1

# SET VARIABLES
SPP=Cbp
REFERENCE=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Cbp2-2_plus_manual_correction.fasta
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/vcf/"$SPP"
TEMP_DIR=/scratch/cfisc004/"$RANDOM"
THREADS=4

mkdir -pv "$TEMP_DIR"

#### PIPELINE #####
# consolidate GVCFs
gatk --java-options "-Xmx16g -Xms4g" GenomicsDBImport \
	--sample-name-map sample_map_"$SPP" \
	--genomicsdb-workspace-path "$RESULTS"/db_"$SLURM_ARRAY_TASK_ID" \
	--tmp-dir "$TEMP_DIR" \
	--reader-threads 3 \
	-XL SCF_1 \
        -XL SCF_2 \
        -XL SCF_3 \
        -XL SCF_4 \
        -XL SCF_5 \
        -XL SCF_6 \
        -XL SCF_7 \
        -XL SCF_8 \
	-XL SCF_9 \
	-XL SCF_10 \
	-XL SCF_11 \
	-XL SCF_12 \
	-XL SCF_13 \
	-XL SCF_14 \
	-XL SCF_15 \
	-XL SCF_16

# joint Genotypying
gatk --java-options "-Xmx16g" GenotypeGVCFs \
   -R "$REFERENCE" \
   -V gendb://"$RESULTS"/db_"$SLURM_ARRAY_TASK_ID" \
   -O "$RESULTS"/"$SLURM_ARRAY_TASK_ID".vcf.gz \
   --tmp-dir "$TEMP_DIR"

# clean up temp
rm -r "$TEMP_DIR"
