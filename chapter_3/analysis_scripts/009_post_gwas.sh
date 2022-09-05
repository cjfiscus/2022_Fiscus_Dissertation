#!/bin/bash

#SBATCH --job-name=post_gwas
#SBATCH -o std/%j.out
#SBATCH -e std/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --ntasks=4
#SBATCH --mem=64gb
#SBATCH -t 1-00:00:00
#SBATCH -p koeniglab

# parse h2
DIR=/rhome/cfisc004/bigdata/projects/cowpea_gwas/results/gwas

rm -f "$DIR"/pve.txt
touch "$DIR"/pve.txt

for i in $DIR/*.log.txt
do
        est=$(grep "pve estimate in the null model" "$i" | cut -f2 -d"=" | xargs)
        nm=$(basename "$i" | sed 's/\.log.txt//g')
        echo -e "$nm""\t""$est" >> "$DIR"/pve.txt
done
##########
# tested snps per gwas
wc -l "$DIR"/*.assoc.txt |  sed 's/\.assoc\.txt//g' | 
	sed 's/\/rhome\/cfisc004\/bigdata\/cowpea_gwas\/results\/gwas\///g' | 
	awk '{print $1,$2}' > "$DIR"/gwas_snps_tested.txt
##########
# export snps with low p
THRES=0.0001
rm -f "$DIR"/gwas_snps_lowp.txt
touch "$DIR"/gwas_snps_lowp.txt

for i in "$DIR"/*.assoc.txt
do
	NAME=$(basename "$i" | sed 's/\.assoc\.txt//g')
	cat "$i" | awk -v THRES="$THRES" -v NAME="$NAME" '$14 < THRES {print NAME,"\t",$0}' >> "$DIR"/gwas_snps_lowp.txt
done
