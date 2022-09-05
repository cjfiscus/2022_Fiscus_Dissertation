#!/bin/bash

DIR=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/gwas/assoc

touch pve.txt
for i in $DIR/*.log.txt
do
	est=$(grep "pve estimate in the null model" "$i" | cut -f2 -d"=" | xargs)
	nm=$(basename "$i" | sed 's/\.log.txt//g')
	echo -e "$nm""\t""$est" >> pve.txt
done
