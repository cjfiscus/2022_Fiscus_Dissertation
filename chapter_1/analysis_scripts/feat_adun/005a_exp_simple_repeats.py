#!/usr/bin/env python3
# Export kmers in list from Kmer table
# cjfiscus
#
# usage:
# python export_kmers.py kmers.txt TABLE
# kmers in kmers.txt must be in first column
# kmers in TABLe must be in second column

import sys
File1=open(sys.argv[1],"r")
File2=open(sys.argv[2],"r")

kmers=[] # kmers to export

# get list of kmers to export
for Line in File1:
	kmer=Line.split()[0]
	kmers.append(kmer)

File1.close()

# export Kmers from table
LineNo = 1 # Line counter
for Line in File2:
	kmer=Line.split()[0]
	
	if kmer in kmers:
		print(Line.strip("\n"))
	else:
		if LineNo == 1:
			print(Line.strip("\n"))
	LineNo = LineNo + 1
