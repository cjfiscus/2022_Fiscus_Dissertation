#!/usr/bin/env Rscript
# filtering

args = commandArgs(trailingOnly=TRUE)

# seq abund table
df<-read.delim(args[1], check.names = F)

# reference 
genes<-read.delim("../../data/A_thaliana_Gene_Database.txt")
repeats<-read.delim("../../data/A_thaliana_Repeat_Database.txt")

# representative non-coding genes
nc_lst<-genes[genes$CLASS=="non coding",]

## noncoding only
noncoding<-df[df$Feature%in%nc_lst$ID,]
nc_rep<-read.delim("../../data/noncoding_rep_loci.txt")

## only representative seqs
noncoding<-noncoding[noncoding$Feature %in% nc_rep$Locus.Identifier,]
##########

# representative pseudogenes
pseudo_lst<-genes[genes$CLASS=="pseudogene",]

## pseudogene only
pseudo<-df[df$Feature%in%pseudo_lst$ID,]
pseudo_rep<-read.delim("../../data/pseudogene_rep_loci.txt")

## only representative seqs
pseudo<-pseudo[pseudo$Feature %in% pseudo_rep$Locus.Identifier,]

##########
## everything else
df1<-df[!df$Feature%in%nc_lst$ID,]
df1<-df1[!df1$Feature%in%pseudo_lst$ID,]

## bind 
df<-as.data.frame(rbind(df1, noncoding, pseudo))
write.table(df, args[2], sep="\t", quote=F, row.names=F)
