#!/usr/bin/env python3
# harmonic mean of beta
# cjfiscus
# 2021-04-13

import sys
import math

file_list=sys.argv[1]

# Function to sum 1/beta for each snp
def sum_beta(file_name):
    File = open(file_name, "r")
    line_num = 1
    for Line in File:
        if line_num > 1:
            parse=Line.split()
            snp_id=parse[0]+":"+parse[2]
            
            beta=1/float(parse[7])
 
            # add to dict
            if snp_id in beta_values.keys():
                beta_old = beta_values.get(snp_id)
                beta_new = beta_old + beta
                beta_values[snp_id] = beta_new
            else:
                beta_values[snp_id] = beta
        line_num=line_num+1
    File.close()
#####

# parse file list
files_to_read=[]
File=open(file_list,"r")
for Line in File:
    parse=Line.strip("\n")
    files_to_read.append(parse)

#####
# process files in filelist
beta_values={}

for i in files_to_read:
    sum_beta(i)

N = len(files_to_read)

# n/sum(1/beta)
for i in beta_values:
    beta_old=beta_values.get(i)
    beta_new = N / beta_old
    beta_values[i] = beta_new

# write out
print("snp_id\tbeta_mean")
for key,val in beta_values.items():
    print(str(key) + "\t" + str(val))
