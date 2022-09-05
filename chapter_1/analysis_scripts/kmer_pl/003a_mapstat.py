## process samtools mapping stat files
## cjfiscus
##
## Output:
## LIBRARY %MAPPED

import re
import sys
import os

# Infile
FileName = sys.argv[1]

File=open(FileName)
Name=FileName.split("/")[-1].strip("_mapstats.txt")

## desired data is on line 5 formatted like this:
# 8513199 + 0 mapped (27.51% : N/A)

LineNo=1

for Line in File:
    if LineNo==5:
        elements=Line.split()
        print(Name + "\t"+ elements[4].strip("(%"))

    LineNo=LineNo+1
