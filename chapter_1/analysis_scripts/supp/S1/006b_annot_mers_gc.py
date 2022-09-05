#!/usr/bin/env python3
## calculate gc content of kmer and append to beginning of line. 
## cjfiscus

## usage
## python annot_mers_gc_content.py InFile.txt OutFile.txt

import sys
import os

LineNum=1 # first line 

File=open(sys.argv[1], "r")

for Line in File:
    ## reset score
    content=0

    ## tabulate count of Gs and Cs
    nts=list(Line.strip("\n").split("\t")[0]) # split into nucleotides

    for nt in nts:
        if str(nt) == "C":
            content = content + 1
        
        elif str(nt) == "G":
            content = content + 1 
            
    ## calculate GC content 
    score = int((content/len(nts))*100)

    ## print score    
    if LineNum > 1:
    	print(str(score)+"\t"+str(Line.strip("\n")))
    
    else: 
    	print("GC_content\t" + str(Line.strip("\n")))
    
    LineNum=LineNum + 1 

