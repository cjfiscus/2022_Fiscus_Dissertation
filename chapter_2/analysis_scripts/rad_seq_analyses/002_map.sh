#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=std/map%j.stdout
#SBATCH --error=std/map%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="map"
#SBATCH -p koeniglab
#SBATCH --array=2-1475

# software dependencies
## bwa-mem2 v.2.0
## samtools 1.10

# load required modules (slurm) 
module load samtools/1.10 

# SET VARIABLES
REFERENCE=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Cbp2-2_plus_manual_correction.fasta
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/bams
SEQLIST=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/sample_map.txt
TEMP_DIR=/rhome/cfisc004/bigdata/temp
THREADS=3

#### PIPELINE #####
# parse sample info
ID=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f1)
SM=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f2)
FILE=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f3)

# mk temp
TEMP_DIR="$TEMP_DIR"/"$ID"
mkdir -pv "$TEMP_DIR"
cd "$TEMP_DIR"

# map
echo "map with bwa-mem2"
bwa-mem2 mem -t "$THREADS" -M "$REFERENCE" \
	-R "@RG\tID:""$ID""\tSM:""$SM""\tPL:ILLUMINA\tLB:""$ID" \
	"$FILE" > "$ID".sam

##########
# sam to sorted bam
echo "samtools sam to bam"
samtools view -b "$ID".sam | samtools sort -T temp - -o "$RESULTS"/"$ID".bam

# clean up temp
rm -r $TEMP_DIR
