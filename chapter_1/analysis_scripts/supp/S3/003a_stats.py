#!/usr/bin/env python3
# analysis.py file1.txt.gz file2.txt.gz
# This script counts the number of total kmers, unique kmers, and singletons. 
# cjfiscus

import sys
import os
import gzip

# process File 1 
File1=gzip.open(sys.argv[1], 'rt')
singletons=0
kmers1=[]

for Line in File1:
    kmer=Line.split("\t")[0]
    count=Line.split("\t")[1].strip("\n")

    kmers1.append(kmer)

    # check if singleton
    if int(count) == 1:
        singletons = singletons + 1

File1.close()

# write out results for file 1
file_name = os.path.split(sys.argv[1])[1].split(".")[0]
Name = "_".join(file_name.split("_")[:-1])
K = file_name.split("_")[-1]

print(Name + "\t" + K + "\t" + "Singletons" + "\t" + str(singletons))
print(Name + "\t" + K + "\t" + "Total" + "\t" + str(len(kmers1)))

# process File 2
File1=gzip.open(sys.argv[2], 'rt')
singletons=0
kmers2=[]

for Line in File1:
    kmer=Line.split("\t")[0]
    count=Line.split("\t")[1].strip("\n")

    kmers2.append(kmer)

    # check if singleton
    if int(count) == 1:
        singletons = singletons + 1

File1.close()

# write out results for file 2
file_name = os.path.split(sys.argv[2])[1].split(".")[0]
Name = "_".join(file_name.split("_")[:-1])
K = file_name.split("_")[-1]

print(Name + "\t" + K + "\t" + "Singletons" + "\t" + str(singletons))
print(Name + "\t" + K + "\t" + "Total" + "\t" + str(len(kmers2)))

# calculate unique kmers
# unique in first set
dif1 = len(list(set(kmers1) - set(kmers2)))
file_name = os.path.split(sys.argv[1])[1].split(".")[0]
Name = "_".join(file_name.split("_")[:-1])
print(Name + "\t" + K + "\t" + "Unique" + "\t" + str(dif1))

# unique in second set
dif2 = len(list(set(kmers2) - set(kmers1)))
file_name = os.path.split(sys.argv[2])[1].split(".")[0]
Name = "_".join(file_name.split("_")[:-1])
print(Name + "\t" + K + "\t" + "Unique" + "\t" + str(dif2))
