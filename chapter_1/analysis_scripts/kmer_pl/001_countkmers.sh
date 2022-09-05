#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --output=std/pl%j.stdout
#SBATCH --error=std/pl%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="redo"
#SBATCH -p koeniglab
#SBATCH --array=924,1134

# software dependencies
## samtools 1.9; 
## trimmomatic 0.36; 
## bedtools 2.27.1; 
## jellyfish 2.2.9; 
## bwa 0.7.17
## mosdepth 0.2.5
## spades 3.15.0

conda deactivate

# load required modules (slurm) 
module load samtools/1.9 trimmomatic/0.36 jellyfish/2.2.9 bwa/0.7.17 

# set python environment (for mosdepth)
module unload miniconda3
module unload miniconda2
module load anaconda3
source activate PyEnv3

# SET VARIABLES
REFERENCE=/rhome/cfisc004/bigdata/projects/K-mers_Arabidopsis/data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
ORGANELLAR=/rhome/cfisc004/bigdata/projects/K-mers_Arabidopsis/data/reference/Arabidopsis_thaliana.TAIR10.Mt.Pt.fa
ADAPTERSPE=/rhome/cfisc004/software/Trimmomatic-0.36/adapters/PE_all.fa
ADAPTERSSE=/rhome/cfisc004/software/Trimmomatic-0.36/adapters/SE_all.fa
DLFAILEDLOG=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/logs/failed.txt
RESULTSDIR=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results
SEQLIST=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/data/1001_Genomes_2016.txt
TEMP_DIR=/scratch/cfisc004
THREADS=8

#### PIPELINE #####

# get filenames from list 
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f15)

# check if sequencing is single end or paired end
if [[ "$FILE" == *";"* ]] ; then
        LIBTYPE="PE" # paired end
else
        LIBTYPE="SE" # single end
fi

# obtain sample name
NAME=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f1)

# make temp directory and go there
TEMP_DIR="$TEMP_DIR"/"$NAME"
mkdir -pv "$TEMP_DIR"
echo "$TEMP_DIR"
cd "$TEMP_DIR"

if [ $LIBTYPE == "PE" ]
then # paired end 
        # Download files
	INDEX=1 # 1 is forward, 2 is reverse
	for i in $(echo $FILE | tr ";" "\n")
	do	
		echo "downloading" "$i" 
		wget -O "$NAME"_"$INDEX".fastq.gz "$i"
		INDEX=$((INDEX + 1))
	done

	# check MD5sums
	echo "verifying checksums..."
	INDEX=1
	MD5SUMS=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f14)
	for i in $(echo "$MD5SUMS" | tr ";" "\n")
	do 
		echo "$i" "$NAME"_"$INDEX".fastq.gz >> chk.md5
    		INDEX=$((INDEX + 1))
	done	

	if md5sum --status -c chk.md5; then
		# continue		
		echo "SUMS GOOD"
	else
		# stop script
		echo "SUMS BAD"
		echo "$NAME" >> "$DLFAILEDLOG"
		rm -r $TEMP_DIR	
		exit 1
	fi

	# Quality/Adapter trimming
	echo "trimming with trimmomatic..."
	java -jar $TRIMMOMATIC PE -threads "$THREADS" \
	"$NAME"_1.fastq.gz "$NAME"_2.fastq.gz \
	"$NAME"_1_trimmed_paired.fq.gz "$NAME"_1_unpaired.fq.gz \
	"$NAME"_2_trimmed_paired.fq.gz "$NAME"_2_unpaired.fq.gz \
	ILLUMINACLIP:"$ADAPTERSPE":2:30:10 \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 HEADCROP:20 MINLEN:36

	# read error correction
	echo "error correction with bayeshammer..."
	spades.py -t "$THREADS" -m 64 --only-error-correction -1 "$NAME"_1_trimmed_paired.fq.gz -2 "$NAME"_2_trimmed_paired.fq.gz -o ./
	mv spades.log "$RESULTSDIR"/logs/error_correct/"$NAME".log

	# map to reference genome
	echo "mapping with bwa..."
	bwa mem -t "$THREADS" -M $REFERENCE ./corrected/"$NAME"_1_trimmed_paired.fq.00.0_0.cor.fastq.gz \
		./corrected/"$NAME"_2_trimmed_paired.fq.00.0_0.cor.fastq.gz > "$NAME"_gen.sam

	# map to organellar genome 
	bwa mem -t "$THREADS" -M $ORGANELLAR ./corrected/"$NAME"_1_trimmed_paired.fq.00.0_0.cor.fastq.gz \
		./corrected/"$NAME"_2_trimmed_paired.fq.00.0_0.cor.fastq.gz > "$NAME"_org.sam

else # single end 
	# Download file
	echo "downloading" "$FILE"	
	wget -O "$NAME".fastq.gz "$FILE"

	# check MD5sum
	echo "verifying checksums..."
	MD5SUMS=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f14)
	echo "$MD5SUMS" "$NAME".fastq.gz >> chk.md5
	
	if md5sum --status -c chk.md5; then
        	# continue
        	echo "SUMS GOOD"
        else
        	# stop script
        	echo "SUMS BAD"
		echo "$NAME" >> "$DLFAILEDLOG"
        	rm -r $TEMP_DIR
		exit 1
        fi

	# Quality/Adapter trimming
	echo "trimming with trimmomatic..."
	java -jar $TRIMMOMATIC SE -threads "$THREADS" \
	"$NAME".fastq.gz "$NAME"_trimmed.fq.gz \
	ILLUMINACLIP:"$ADAPTERSSE":2:30:10 \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 HEADCROP:20 MINLEN:36 

	# read error correction
	echo "error correction with bayeshammer..."
	spades.py -t "$THREADS" -m 64 --only-error-correction -s "$NAME"_trimmed.fq.gz -o ./
	mv spades.log "$RESULTSDIR"/logs/error_correct/"$NAME".log

	# map to reference genome
	echo "mapping with bwa..."
	bwa mem -t "$THREADS" -M $REFERENCE ./corrected/"$NAME"_trimmed.fq.00.0_0.cor.fastq.gz  > "$NAME"_gen.sam
	
	# map to organellar genome
	bwa mem -t "$THREADS" -M $ORGANELLAR ./corrected/"$NAME"_trimmed.fq.00.0_0.cor.fastq.gz  > "$NAME"_org.sam

fi

# sam to sorted bam
echo "samtools sam to bam"
# only keep primary alignments for ref genome mapping (-F 256)
samtools view -F 256 -bS "$NAME"_gen.sam | samtools sort -T temp_Pt - -o "$NAME"_gen.bam

# only keep reads that did not map to organelles
samtools view -f4 -bS "$NAME"_org.sam | samtools sort -T temp_Pt - -o "$NAME".unmapped.bam

# mapping stats
echo "samtools flagstat"
samtools flagstat "$NAME"_gen.bam > $RESULTSDIR/mappingstats/genome/"$NAME"_mapstats.txt
samtools flagstat "$NAME"_org.sam > $RESULTSDIR/mappingstats/organellar/"$NAME"_mapstats.txt

# index bam
#echo "samtools indexing bam"
samtools index "$NAME"_gen.bam

# calculate coverage of ref per base and in 1kb windows
echo "calculating coverage with mosdepth..."
mosdepth -t "$THREADS" -b 1000 "$RESULTSDIR"/coverage/genome/"$NAME" "$NAME"_gen.bam

# subset reads that did not map to organelles
echo "extracting unmapped reads..."
# export unmapped reads
bedtools bamtofastq -i "$NAME".unmapped.bam -fq "$NAME".unmapped.fq

# Count 12-mers in reads that did not map to organelles
echo "counting K-mers with jellyfish"
jellyfish count -C -m 12 -s 500M -t "$THREADS" -o "$NAME".jf "$NAME".unmapped.fq
jellyfish dump -tc "$NAME".jf > $RESULTSDIR/kmer_counts/"$NAME".txt

# clean up
rm -r $TEMP_DIR
