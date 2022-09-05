#!/usr/bin/env Rscript

# script to normalize seq. abundance with buscos

# read sysargs
args = commandArgs(trailingOnly=TRUE)
# args[1] is abund table
# args[2] is busco table
# args[3] is norm abund table
# args[4] is norm busco table 

# load pkg 
library(pacman)
p_load(matrixStats)
options(stringsAsFactors = F)

# read in seq. abund ests. 
abund<-read.delim(args[1], check.names=F)

## transform
abund[2:ncol(abund)]<-2^abund[2:ncol(abund)]

# read in busco abund ests. 
busco<-read.delim(args[2], check.names=F)

## transform
busco[2:ncol(busco)]<-2^busco[2:ncol(busco)]

# calculate median abund. est. for buscos within each lib
meds<-as.data.frame(cbind(names(busco)[2:ncol(busco)],
                          colMedians(as.matrix(busco[2:ncol(busco)]))))
names(meds)<-c("acc", "med")
meds$med<-as.numeric(as.character(meds$med))

# normalize abund. ests. for buscos
for (i in 2:ncol(busco)){
  index<-i-1
  busco[,i]<-busco[,i]/meds[index,2]
}

# normalize abund. ests. for all other sequences
for (i in 2:ncol(abund)){
  index<-i-1
  abund[,i]<-abund[,i]/meds[index,2]
}

## write out tables
abund1<-as.data.frame(cbind(abund[,1], round(abund[,2:ncol(abund)],4)))
names(abund1)[1]<-"Feature"
write.table(abund1, file=args[3], sep="\t", quote=F, row.names=F)

busco1<-as.data.frame(cbind(busco[,1], round(busco[,2:ncol(busco)],4)))
names(busco1)[1]<-"Feature"
write.table(format(busco1, digits=4, scientific=F), file=args[4], sep="\t", quote=F, row.names=F)
