#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=200gb
#SBATCH --output=pp%j.stdout
#SBATCH --error=pp%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="post"
#SBATCH -p koeniglab

## TEMP is dir in RAM, WKDIR is directory containing individual K-mer count files, OUT is final merged count file

TEMP=/dev/shm
WKDIR=../results/counts
OUT=/rhome/cfisc004/bigdata/projects/SX/results/table_12mers.txt

## set working directory (where ind. K-mer count files are) 
cd $WKDIR

## join the first two files to make output 
FILE1="$(ls | head -n 1)" # first file
FILE2="$(ls | head -n 2 | tail -n 1)" # second file 

## make a list of the remaining files to iterate through
ls *.txt.gz | tail -n +3 > $TEMP/list."$SLURM_JOB_ID".tmp # does not include 1st 2 files

## sort file by K-mer 
sort <(zcat "$FILE1") > $TEMP/"$FILE1.tmp"
sort <(zcat "$FILE2") > $TEMP/"$FILE2.tmp"

## merge first two files to create master
join -1 1 -2 1 -a 1 -a 2 -e 0 -o 0,1.2,2.2 -t $'\t' $TEMP/"$FILE1.tmp" $TEMP/"$FILE2.tmp" > $TEMP/"$SLURM_JOB_ID".txt

## write header line 
NAME1=$(basename "$FILE1" | cut -d. -f1) # strip extension from filenames 
NAME2=$(basename "$FILE2" | cut -d. -f1)
HEADER=$(echo "mer""\t""$NAME1""\t""$NAME2")

## remove temporary files
rm $TEMP/"$FILE1.tmp"; rm $TEMP/"$FILE2.tmp" 

NUM="$(ls *.txt.gz | tail -n +3 | wc -l)" # number of remaining files 

for i in `seq 1 $NUM` # do this for remaining files 
	do 
		FILE=$(head -n $i $TEMP/list."$SLURM_JOB_ID".tmp | tail -n 1)
		
		## sort file by K-mer 
		sort <(zcat "$FILE") > $TEMP/"$FILE.tmp" 
		
		## Merge with master table
		join -1 1 -2 1 -a 1 -a 2 -e 0 -o auto -t $'\t' $TEMP/"$SLURM_JOB_ID".txt $TEMP/"$FILE.tmp" > $TEMP/temp."$SLURM_JOB_ID".txt
		mv $TEMP/temp."$SLURM_JOB_ID".txt $TEMP/"$SLURM_JOB_ID".txt 
		
		## append Filename to header line 
		NAME=$(basename "$FILE" | cut -d. -f1) # grab filename
		HEADER="$(echo "$HEADER""\t""$NAME")"
		
		## remove temporary files
		rm $TEMP/"$FILE.tmp"

	done 

## add header line to file 
sed '1s/^/'"$HEADER"'\n/' $TEMP/"$SLURM_JOB_ID".txt > $TEMP/temp."$SLURM_JOB_ID".txt
mv $TEMP/temp."$SLURM_JOB_ID".txt $OUT 

## remove temporary files 
rm $TEMP/list."$SLURM_JOB_ID".tmp $TEMP/"$SLURM_JOB_ID".txt
