#!/usr/bin/env Rscript
# initial kmer table filtering
# cjfiscus

library(pacman)
p_load(data.table)

args = commandArgs(trailingOnly=TRUE)

# read in data 
counts<-fread(args[1])

# id samples to be filtered
drop<-read.table("../results/filters/filtered1.txt")
length(drop)

# list of colnames to retain
keep<-names(counts)[!(names(counts) %in% drop$V1)]

# filter table
counts<-counts[ , ..keep]
dim(counts)

fwrite(counts, file=args[2], sep="\t", quote=F)
