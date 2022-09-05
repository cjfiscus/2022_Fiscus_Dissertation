#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=400gb
#SBATCH --output=pp%j.stdout
#SBATCH --error=pp%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="post2"
#SBATCH -p koeniglab

## TEMP is dir in RAM, WKDIR is directory containing individual K-mer count files, OUT is final merged count file

TEMP=/dev/shm
WKDIR=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/kmer_counts2
OUT=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/athaliana_12mers2.txt

## set working directory (where ind. K-mer count files are) 
cd $WKDIR

# define function for multiple join
joinm() {
f1=$1; f2=$2; shift 2;
if [ $# -gt 0 ]
then
# join files
join -a 1 -a 2 -e 0 "$f1" "$f2" -o auto -t $'\t' | joinm - "$@"

else
# join 2 files
join -a 1 -a 2 -e 0 "$f1" "$f2" -o auto -t $'\t'

fi
}

## merge files
joinm *.txt > "$OUT".tmp

## construct header
HEADER=$(echo "mer")

NUM=$(ls *.txt | wc -l)
for i in `seq 1 $NUM`
do
NAME=$(ls *.txt | head -n $i | tail -n1)
NAME=$(basename "$NAME" | cut -d. -f1)
HEADER=$(echo "$HEADER""\t""$NAME")
done

# add header to file
sed '1s/^/'"$HEADER"'\n/' "$OUT".tmp > "$OUT"

## cleanup
rm "$OUT".tmp
