## GC Content Calculator for K-mer data
## cjfiscus
#
# Input is K-mer count table where K-mer is in first col
# Output is list of Library\tGC_content
#
## Usage: 
## 004d_calc_gc_lib.py InFile.txt OutFile.txt 

import sys
import os

# Infile
File = sys.argv[1]
File1 = sys.argv[2]

InFile=open(File,"r")
OutFile=open(File1,"w")

names=[] # ({index:name})
GC_sum=[]
total_sum=[]

LineNum=1

## build a list of names 
for Line in InFile:
    if LineNum == 1:
        elements=Line.strip("\n")
        elements=elements.split("\t")
        
        for i in range(1,len(elements)):
            names.append(elements[i])
            GC_sum.append(0) ## initialize lists 
            total_sum.append(0)

    else:
        pass  

    LineNum=LineNum + 1
InFile.close()

## take counts 
InFile=open(File,"r")
LineNum = 1

for Line in InFile:
    if LineNum > 1:
        elements=Line.strip("\n")
        elements=elements.split("\t")

        mer=list(elements[0])
        ## GC sum
        if mer[0] == "G" or mer[0] == "C":
            for i in range(len(GC_sum)):
                GC_sum[i]= float(GC_sum[i]) + float(elements[1+i])
                print(GC_sum[i])
        else:
            pass 
        
        ## total sum     
        for i in range(len(total_sum)):
            total_sum[i] =float(total_sum[i]) + float(elements[1+i]) 

    else: # skip first line 
        pass 

    LineNum=LineNum + 1

InFile.close()

## Write results to OutFile
OutFile.write("Library\tGC_content\n")

for i in range(0,len(GC_sum)):
    
    gc_content=(GC_sum[i]/total_sum[i])*100
    gc_content="{0:.2f}".format(gc_content)

    OutFile.write(str(names[i])+"\t"+str(gc_content)+"\n")
