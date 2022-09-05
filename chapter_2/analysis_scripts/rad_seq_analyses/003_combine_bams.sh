#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00
#SBATCH --job-name="map"
#SBATCH -p koeniglab
#SBATCH --array=2-1475

## ARRAY 2-1475

# software dependencies
## samtools 1.10

# load required modules (slurm) 
module load samtools/1.10 

FILE=/rhome/cfisc004/bigdata/projects/capsella_radseq/data/sample_map.txt
SM=$(head -n "$SLURM_ARRAY_TASK_ID" "$FILE" | tail -n1 | cut -f2)
TMP="$SLURM_ARRAY_TASK_ID"_tmp.txt
OUT=/rhome/cfisc004/bigdata/projects/capsella_radseq/results/processed_bams/"$SM".bam

## pull lst of bams
grep "$SM" "$FILE" > "$TMP"

IN=$(cut -f1 "$TMP" | sed 's/^/\/rhome\/cfisc004\/bigdata\/projects\/capsella_radseq\/results\/bams\//g' | sed 's/$/\.bam/g')

# check if 1 bam or multiple bams
if [ $(wc -l "$TMP" | cut -d" " -f1 | bc) -gt 1 ] 
then # merge multiple bams
readarray -t ARR <<< "$IN"
echo "merging ${ARR[@]} into $OUT"
samtools merge "$OUT" "${ARR[@]}"

else # rename single bam
echo "cp $IN $OUT"
cp "$IN" "$OUT"
fi

# cleanup tmp file
rm $TMP
