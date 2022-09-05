#!/usr/bin/env python3
# faidx_sum.py
# cjfiscus
# 7/16/2018
#
# This script takes a samtools faidx index and returns the sum of seq lengths.  

import sys
import os

File=open(sys.argv[1], "r")

seqlen = 0

for Line in File:
    add = Line.split()[1]
    seqlen = int(seqlen) + int(add)

FileName=os.path.split(sys.argv[1])[1].split(".")[0]
print(str(FileName) + "\t" + str(seqlen))

File.close()
