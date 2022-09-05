#!/usr/bin/env python3
# Non-overlapping windows from .fai
# cjfiscus
# 9/13/19
#
# usage:
# python gen_intervals.py genome.fai 10000 

import sys

# parse args
# samtools faidx index
File = open(sys.argv[1], "r")

# size of interval
size = int(sys.argv[2])

for Line in File:
    # get chr and length for chr
    parse = Line.strip("\n").split()
    chrom = parse[0]
    length = parse[1]
    
    # how many windows to generate
    num_windows = int(length)//int(size)

    # chr shorter than window
    if num_windows == 0:
        print(str(chrom) + ":1-" + str(length))

    else: # print ranges 
        for i in range(1, num_windows + 1):
            start = i * size - size + 1
            end = i * size
            print(str(chrom) + ":" + str(start) +"-" + str(end))
