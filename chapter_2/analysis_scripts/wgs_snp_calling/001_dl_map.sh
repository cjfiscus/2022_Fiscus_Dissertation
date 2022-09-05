#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=32G
#SBATCH --output=std/map%j.stdout
#SBATCH --error=std/map%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --job-name="co"
#SBATCH -p koeniglab
#SBATCH --array=2-4

# software dependencies
## bwa-mem2 v.2.0
## samtools 1.10

# load required modules (slurm) 
module load samtools/1.10 

# SET VARIABLES
REFERENCE=/rhome/cfisc004/shared/GENOMES/CAPSELLA/FINISHED_v2/data/assemblies/Co39_plus_manual_correction.fasta
RESULTS=/rhome/cfisc004/bigdata/projects/capsella_genomes/results/alignments/Co
SEQLIST=../../data/Capsella_WGS_Sample_List.txt
TEMP_DIR=/scratch/cfisc004
THREADS=4

#### PIPELINE #####
# parse sample file
NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f4)
FILES=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1| cut -f6)
MD5SUMS=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f5)
ID=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f2)

# mk temp
TEMP_DIR="$TEMP_DIR"/"$NAME"_"$RANDOM"
mkdir -pv "$TEMP_DIR"
cd "$TEMP_DIR"

# check if sequencing is single end or paired end
if [[ "$FILES" == *";"* ]] ; then
        LIBTYPE="PE" # paired end
else
        LIBTYPE="SE" # single end
fi

# paired end
if [ $LIBTYPE == "PE" ] 
then 
	# dl fqs
	INDEX=1
	for i in $(echo "$FILES" | tr ";" "\n")
	do 
		echo "downloading" "$i"
		wget -O "$NAME"_"$INDEX".fq.gz "$i"
		INDEX=$((INDEX + 1))
	done

	# check md5sums
	echo "verifiying checksums..."
	INDEX=1
	for i in $(echo "$MD5SUMS" | tr ";" "\n")
	do
		echo "$i" "$NAME"_"$INDEX".fq.gz >> chk.md5
		INDEX=$((INDEX + 1))
	done

	if md5sum --status -c chk.md5; then
		echo "SUMS GOOD"
	else
		echo "SUMS BAD"
		exit 1
	fi 

	# map
	echo "map with bwa-mem2"
	bwa-mem2 mem -t "$THREADS" -M "$REFERENCE" \
		-R "@RG\tID:""$ID""\tSM:""$NAME""\tPL:ILLUMINA\tLB:""$ID" \
		"$NAME"_1.fq.gz "$NAME"_2.fq.gz > "$NAME".sam


else # single end
	# dl fqs
	echo "downloading" "$FILES"
	wget -O "$NAME".fq.gz "$FILES"

	# check md5sums
	echo "verifiying checksums..."
	echo "$MD5SUMS" "$NAME".fq.gz >> chk.md5

	if md5sum --status -c chk.md5; then
		echo "SUMS GOOD"
	else
		echo "SUMS BAD"
		exit 1
	fi 

	# map
	echo "map with bwa-mem2"
	bwa-mem2 mem -t "$THREADS" -M "$REFERENCE" \
		-R "@RG\tID:""$ID""\tSM:""$NAME""\tPL:ILLUMINA\tLB:""$ID" \
		"$NAME".fq.gz > "$NAME".sam
fi

##########
# sam to sorted bam
echo "samtools sam to bam"
## this line for ref only
#samtools view -b "$NAME".sam | samtools sort -T temp - -o "$RESULTS"/"$ID".bam

## this line for everything else
samtools view -b "$NAME".sam | samtools sort -T temp - -o "$RESULTS"/"$NAME".bam

# clean up temp
rm -r $TEMP_DIR
