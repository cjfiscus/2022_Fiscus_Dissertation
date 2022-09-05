#!/usr/bin/env Rscript
# pca with rsvd
# cjfiscus

library(pacman)
p_load(data.table, rsvd)

args = commandArgs(trailingOnly=TRUE)

# read in GC corrected, normalized data
print("reading in data")
counts<-fread(args[1])
dim(counts)
counts[1:5,1:5]

## DROP OUTGROUPS
drop<-c("A_lyrata1", "A_lyrata2", "C_rubella1", "C_rubella2")
counts[,c(drop):=NULL]
dim(counts)

# prepare data for PCA 
print("formatting table")
counts<-as.data.frame(counts)
counts[is.na(counts)]<-0
row.names(counts)<-counts$mer 
counts$mer<-NULL

# determine minimum non-zero value (to add to 0s for log)
minValue<-min(counts[counts > 0])

# log transform
counts<-log2(counts+minValue)

# transpose (samples will be rownames)
Norm.t<-t(counts)
rm(counts) # purge from RAM 

# pca
print("pca")
df.pca<-rpca(Norm.t, k=20, center=T, scale=F)

summary(df.pca) 

# percent variance for each component
per.var<-as.data.frame(cbind(1:10,summary(df.pca)[3,1:10])) # calculate % variance 
colnames(per.var)<-c("PC", "PercentVariance")
out<-paste0("../results/pca/pervar_", args[2], ".txt.gz")
write.table(per.var, gzfile(out), sep="\t", quote=F,row.names=F) # write table out 

## get coords out
coord<-as.data.frame(df.pca$x) ## coordinates
out<-paste0("../results/pca/coords_", args[2], ".txt.gz")
write.table(coord, gzfile(out), quote=F, sep="\t")

## get rotations out 
rot<-as.data.frame(df.pca$rotation)
out<-paste0("../results/pca/rot_", args[2], ".txt.gz")
write.table(rot, gzfile(out), quote=F, sep="\t")

## get eigenvalues out
ev<-as.data.frame(df.pca$eigvals)
out<-paste0("../results/pca/eigvals_", args[2], ".txt.gz")
write.table(ev, gzfile(out), quote=F, sep="\t")

# get sdev out
sdev<-as.data.frame(df.pca$sdev)
write.table(sdev, gzfile("../results/pca/sdev.txt.gz"), quote=F, sep="\t")
