#!/usr/bin/env python3
# Calculate X^2 using Fisher's Method
# cjfiscus
# 11/20/19
# modified 03/23/2021

import sys
import math

file_list=sys.argv[1]

# Function to sum ln(p) for each snp
def sum_chi(file_name):
    File = open(file_name, "r")
    line_num = 1
    for Line in File:
        if line_num > 1:
            parse=Line.split()
            snp_id=parse[0]+":"+parse[2]
            
            # check if p value can be represented
            if float(parse[13]) > 0:
                chi=math.log(float(parse[13]))
            
            # if too low just set to low p value
            else:
                chi=math.log(1e-100)

            # add to dict
            if snp_id in chi_values.keys():
                chi_old = chi_values.get(snp_id)
                chi_new = chi_old + chi
                chi_values[snp_id] = chi_new
            else:
                chi_values[snp_id] = chi
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
chi_values={}

for i in files_to_read:
    sum_chi(i)

# multiply by -2
for i in chi_values:
   chi_old = chi_values.get(i)
   chi_new = -2 * chi_old
   chi_values[i] = chi_new

# write out
print("snp_id\tchisq")
for key,val in chi_values.items():
    print(str(key) + "\t" + str(val))
