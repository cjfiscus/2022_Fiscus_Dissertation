#!/usr/bin/env Rscript
# calc basic stats
# cjfiscus
# 2022-03-29

library(pacman)
p_load(data.table, readxl, adegenet, hierfstat, pegas, stringr)

# import genotypes and filter
df1<-read_excel("../data/IITA_Core_FinalReport.xlsx", sheet="2021 Core Accessions")
df1<-df1[4:nrow(df1),]
df1<-as.data.frame(df1)
names(df1)<-df1[1,]
df1<-df1[2:nrow(df1),]
df1<-df1[complete.cases(df1),]

## rm sites not in nuclear genome
df1<-df1[!df1$`chrom 1.0`=="XX",]

## subset to kept sites
sites<-read.table("../data/cowpea_gwas.bim")
df1<-df1[df1$SNP_Name %in% sites$V2,]
df1<-df1[,-c(2,3)]

## subset to kept samples
row.names(df1)<-df1$SNP_Name
df1$SNP_Name<-NULL
names(df1)<-str_replace(names(df1), "Tvu","TVu")
samps<-read.table("../data/samples.txt")
df1<-df1[,names(df1) %in% samps$V1]

# convert to obj
df1<-as.data.frame(t(df1))
ind<-as.character(row.names(df1))
population<-as.character(rep(1,nrow(df1)))
geno<-df2genind(df1, ploidy=2, ind.names=ind, pop=population, sep="", ncode=1,
		NA.char="-")

# calc Ho and He
div<-summary(geno)
het<-as.data.frame(cbind(div$Hexp, div$Hobs))
names(het)<-c("He", "Ho")
het$snp_id<-row.names(het)
het<-het[,c("snp_id", "He", "Ho")]
write.table(het, "../results/stats/het.txt", sep="\t", quote=F, row.names=F)

## test for diff b/t He and ho
test<-bartlett.test(list(div$Hexp, div$Hobs))
htest<-as.data.frame(cbind(test$statistic, test$parameter, test$p.value))
names(htest)<-c("K2", "df", "p-value")
write.table(htest, "../results/stats/Ho_He_bartlett.txt", sep="\t", quote=F, row.names=F)
##########

# calc af
system("plink --bfile ../data/cowpea_gwas --freq --allow-extra-chr --out ../results/stats/af")
system("plink --bfile ../data/cowpea_gwas --within ../data/strata_bio.txt --freq --allow-extra-chr --out ../results/stats/af")

# calc pw dist
system("plink --bfile ../data/cowpea_gwas --allow-extra-chr --distance square --out ../results/stats/dist")

# calc F
## prune and filter snps
system("plink --bfile ../data/cowpea_gwas --indep-pairwise 50 10 0.2 --maf 0.01 --geno 0.05 --allow-extra-chr --out ../results/stats/prune")
system("plink --bfile ../data/cowpea_gwas --extract ../results/stats/prune.prune.in --allow-extra-chr 0 --make-bed --out ../results/stats/pruned")
system("plink --bfile ../results/stats/pruned --allow-extra-chr --het --out ../results/stats/fstat")
