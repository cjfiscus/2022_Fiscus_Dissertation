#!/usr/bin/env Rscript
# this script produces phenotype files ready for gwas
# cjfiscus

options(stringsAsFactors = FALSE)
library(pacman)
p_load(stringr, dplyr, tidyr)

# import accession list
samples<-read.table("../../results/gwas/id_lst.txt")
names(samples)<-"id"

# import sequence list
seqs<-read.table("../../data/variable_sequences.txt")
names(seqs)<-"Feature"

# import PC coordinates
coords<-read.delim(gzfile("../../results/pca/coords_qq_corrected.txt.gz"))
coords<-coords[,1:3]
names(coords)<-c("PC1_kmer", "PC2_kmer", "PC3_kmer")
coords<-as.data.frame(cbind(row.names(coords), coords))
names(coords)[1]<-"id"

# import seq abundance estimates
seq_abund<-read.delim("../../results/feat_abund/A_thal_kmer_abund_norm.txt")

# subset variable seqs only 
seq_abund<-seq_abund[seq_abund$Feature %in% seqs$Feature,]
rownames(seq_abund)<-seq_abund$Feature
seq_abund$Feature<-NULL

# transpose table
seq_abund<-as.data.frame(t(seq_abund))
seq_abund$id<-row.names(seq_abund)
seq_abund$id<-str_remove(seq_abund$id, "X")

# import and format dist
dist1<-read.table("../../results/dist/dist_anc_noscale.txt")
names(dist1)<-c("id", "anc", "value")
dist1<-spread(dist1, anc, value)


dist1$dist_lyrata<-(dist1$A_lyrata1+dist1$A_lyrata2)/2
dist1$dist_rubella<-(dist1$C_rubella1+dist1$C_rubella2)/2
dist1<-dist1[,c("id", "dist_lyrata", "dist_rubella")]

# 
# join tables 
pheno<-merge(coords, dist1, by="id")
pheno<-merge(pheno, seq_abund, by="id")

# subset accessions
pheno<-pheno[pheno$id %in% samples$id,] 

# find minimum value (for log transformation)
minValue<-min(pheno[,7:ncol(pheno)][pheno[,7:ncol(pheno)] > 0])

# transform sequence abundance estimates (leave PC coords and dist alone) 
pheno[,7:ncol(pheno)]<-log2(pheno[,7:ncol(pheno)] + minValue)

pheno<-pheno[order(as.numeric(as.character(pheno$id))),]

# write out phenotypes
# write out transformed phenotypes 
write.table(pheno, "../../results/gwas/phenotypes_transformed.txt", sep="\t", quote=F, col.names=T, row.names=F)
