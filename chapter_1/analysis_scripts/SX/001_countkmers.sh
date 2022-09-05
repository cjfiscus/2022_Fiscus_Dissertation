#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=pl%j.stdout
#SBATCH --error=pl%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="pl"
#SBATCH -p batch
#SBATCH --array=2-5

# software dependencies
## axel 2.16.1 (https://github.com/axel-download-accelerator/axel/releases) 
## jellyfish 2.2.9; 

# load required modules (slurm) 
module load jellyfish/2.2.9

# SET VARIABLES
WORKINGDIR=../
RESULTS_DIR=/rhome/cfisc004/bigdata/projects/SX/results/
SEQLIST=/rhome/cfisc004/bigdata/projects/SX/data/Michael_2018.txt
TEMP_DIR=../data/
THREADS=4

#### PIPELINE #####

# get filenames from list 
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f5)

# check if sequencing is single end or paired end
if [[ "$FILE" == *";"* ]] ; then
        LIBTYPE="PE" # paired end
else
        LIBTYPE="SE" # single end
fi

# obtain sample name
NAME=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f3)

# make temp directory and go there
TEMP_DIR="$TEMP_DIR""$NAME"
mkdir -pv "$TEMP_DIR"
cd "$TEMP_DIR"

if [ $LIBTYPE == "PE" ]
then # paired end 
        # Download files
	INDEX=1 # 1 is forward, 2 is reverse
	for i in $(echo $FILE | tr ";" "\n")
	do	
		echo "downloading" "$i" 
		#axel -n "$THREADS" "$i" -o "$NAME"_"$INDEX".fastq.gz
		INDEX=$((INDEX + 1))
	done

	# check MD5sums
	echo "verifying checksums..."
	INDEX=1
	MD5SUMS=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f4)
	for i in $(echo "$MD5SUMS" | tr ";" "\n")
	do 
		echo "$i" "$NAME"_"$INDEX".fastq.gz >> chk.md5
    		INDEX=$((INDEX + 1))
	done	

	if md5sum --status -c chk.md5; then
		# continue		
		echo "SUMS GOOD"

		# count k-mers 
		jellyfish count -C -m 12 -s 500M -t 4 -o "$RESULTS_DIR"/"$NAME".jf -F 2 <(zcat "$NAME"_1.fastq.gz) <(zcat "$NAME"_2.fastq.gz)
		jellyfish dump -tc "$RESULTS_DIR"/"$NAME".jf | gzip > "$RESULTS_DIR"/"$NAME".txt.gz	
	else
		# stop script
		echo "SUMS BAD"
		rm -r $TEMP_DIR	
		exit 1
	
		
	fi

else # single end 
	# Download file
	echo "downloading" "$FILE"	
	#axel -n "$THREADS" "$FILE" -o "$NAME".fastq.gz
	
	# check MD5sum
	echo "verifying checksums..."
	MD5SUMS=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f4)
	echo "$MD5SUMS" "$NAME".fastq.gz >> chk.md5
	
	if md5sum --status -c chk.md5; then
        	# continue
        	echo "SUMS GOOD"

		# count k-mers
		jellyfish count -C -m 12 -s 500M -t 4 -o "$RESULTS_DIR"/"$NAME".jf <(zcat "$NAME".fastq.gz)
		jellyfish dump -tc "$RESULTS_DIR"/"$NAME".jf | gzip > "$RESULTS_DIR"/"$NAME".txt.gz
        else
        	# stop script
        	echo "SUMS BAD"
        	rm -r $TEMP_DIR
		exit 1
        fi

fi
