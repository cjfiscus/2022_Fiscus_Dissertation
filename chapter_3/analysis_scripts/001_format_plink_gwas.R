#!/usr/bin/env Rscript
# format cowpea genos and phenos into plink text
# cjfiscus
# 2021-11-09
# updated 2022-04-05
# 
# FILTERING APPLIED
# filter sites w/ > 95% missing data
# filter samples w/ > 50% missing data
# filter samples w/ > 20% het calls
# filter near identical lines using Hamming dist threshold of 500
##########

library(pacman)
p_load("readxl", "stringr", "tidyr", "data.table", "fastDummies", "ggplot2", "cowplot", "ggdendro")

# FORMAT GENOTYPES INTO TPED
## read in sheet 1
df1<-read_excel("../data/IITA_Core_FinalReport.xlsx", sheet="2021 Core Accessions")

## shape up data
### remove extra rows at top
df1<-df1[4:nrow(df1),]
df1<-as.data.frame(df1)
names(df1)<-df1[1,]
df1<-df1[2:nrow(df1),]

### remove extra rows at bot
df1<-df1[complete.cases(df1),]

### rm sites not in nuclear genome
df1<-df1[!df1$`chrom 1.0`=="XX",]

# gentle filter
## calc prop missing for sites
missSite<-as.data.frame(cbind(df1$SNP_Name, rowSums(df1=="--")))
missSite$prop<-as.numeric(missSite$V2)/length(names(df1[4:ncol(df1)]))
names(missSite)<-c("snp_id", "n", "prop_miss")
quantile(missSite$prop_miss)
p1<-ggplot(missSite, aes(x=prop_miss)) + geom_density() +
  theme_cowplot() +
  xlab("prop. calls missing per site")
ggsave("prop_miss_site.pdf", p1, height=4, width=4)

write.table(missSite, "missing_per_site.txt", quote=F, row.names=F)

## filter sites w/ > 95% missing data
keep<-missSite[missSite$prop <=0.95,]
df1<-df1[df1$SNP_Name %in% keep$snp_id,]

## calc prop miss per smp
missSamp<-as.data.frame(cbind(names(df1), colSums(df1=="--")))
missSamp<-missSamp[4:nrow(missSamp),]
missSamp$V2<-as.numeric(missSamp$V2)
missSamp$prop<-missSamp$V2/nrow(df1)
quantile(missSamp$prop)
p1<-ggplot(missSamp, aes(x=prop)) + geom_density() + theme_cowplot() +
  xlab("prop. calls missing per sample") 
ggsave("prop_miss_samp.pdf", p1, height=4, width=4)

names(missSamp)<-c("sample", "n", "prop")
write.table(missSamp, "missing_per_sample.txt", sep="\t", quote=F, row.names=F)
##########

## filter samples w/ > 50% missing data
filt<-missSamp[missSamp$prop > 0.50,]$sample
filt
df1<-as.data.frame(df1)
df1<-df1[,!names(df1) %in% filt]

##########
# create tped formatted
## format columns 1-4
info<-df1[,1:3]
info$dummy<-0
info<-info[,c("chrom 1.0", "SNP_Name", "dummy", "Position")]
names(info)<-c("chr", "id", "dummy", "bp")

## format genotypes
calls<-df1[,4:ncol(df1)]
names(calls)<-str_replace(names(calls), "Tvu", "TVu")
samps<-names(calls)

## set missing to 0
calls<-as.data.frame(apply(calls, 2, str_replace, pattern="--", replacement="00"))

## split cols 
df2 <- as.data.frame(unlist(lapply(calls, data.table::tstrsplit, ""),
                              recursive = FALSE))

## combine with first four cols
out<-as.data.frame(cbind(info, df2))

## cleanup
rm(df1); rm(df2)

## write out table
write.table(out, "../data/cowpea.tped", sep=" ", quote=F, row.names=F, col.names=F)
##########

# FORMAT PHENOTYPES INTO TFAM
## phenotypes
## read in data
df1<-read_excel("../data/ID-TVu CrossReference (Tchamba_11August2020).xls", sheet="Cowpea characterization")
df2<-read_excel("../data/IITA-Cowpea collection-Cowpea Evaluation.xls")
df1<-merge(df1, df2, by="ID")
df1<-df1[df1$`Accession name` %in% samps,]
df1$`Thrips screening`<-NULL
df1$`Pod bugs screening`<-NULL

## generate list of columns to make into dummies (encode qual traits)
dums<-as.data.frame(cbind(names(df1), unlist(lapply(df1, is.character))))
dums<-dums[dums$V2=="TRUE",]
df<-df1[,colnames(df1) %in% dums$V1]
dums<-dums[!dums$V1=="Accession name",]
dums<-dums$V1

## subset numeric data & log
df2<-df1[,!colnames(df1) %in% dums]
df2$ID<-NULL
df2[,2:ncol(df2)]<-log(df2[,2:ncol(df2)] + min(df2[,2:ncol(df2)], na.rm = T))

## code vars w/ dummy vars
df<-dummy_cols(df, select_columns=dums, ignore_na=T)
df<-df[,!colnames(df) %in% dums]

## merge data
df<-merge(df2, df, by="Accession name")

## match ids with genotype data
lst<-as.data.frame(cbind(samps, 1))
names(lst)[1]<-"Accession name"
df2<-merge(lst, df, by=c("Accession name"), all.x=T)
df2$V2<-NULL
df2<-df2[match(samps,df2$`Accession name`),]
all.equal(df2$`Accession name`, samps)

# format output
## add first five cols
out2<-as.data.frame(cbind(df2$`Accession name`, df2$`Accession name`, 0, 0, 0))
out2<-as.data.frame(cbind(out2, df2[,2:ncol(df2)]))

## format phenotype names
phenos<-as.data.frame(names(out2)[6:ncol(out2)])
phenos$N<-1:nrow(phenos)
names(phenos)<-c("trait", "N")
phenos<-phenos[,c("N", "trait")]

## remove duplicate phenotypes
### define dups
ind<-c(seq(105,143), seq(198,251))
dups<-phenos[row.names(phenos) %in% ind,2]

### rm dups from output
out2<-out2[,!names(out2) %in% dups]
phenos<-as.data.frame(names(out2)[6:ncol(out2)])
phenos$N<-1:nrow(phenos)
names(phenos)<-c("trait", "N")
phenos<-phenos[,c("N", "trait")]

# write out data
write.table(out2, "../data/cowpea.tfam", sep=" ", quote=F, row.names=F, col.names=F)
write.table(phenos, "../data/phenotypes.txt", 
            quote=F, row.names=F, col.names=F, sep="\t")

##########
## calc het per sample for filtering
system("plink --tfile ../data/cowpea --allow-extra-chr --het")

### import data
het<-read.table("plink.het", header=T)
het$prop_het<-(het$`N.NM.` - het$`O.HOM.`)/(het$`N.NM.`)
quantile(het$prop_het)
p1<-ggplot(het, aes(x=prop_het)) + geom_density() +
	theme_cowplot()
ggsave("het_per_samp.pdf", p1, height=4, width=4)

### lst of samps > 0.20 het calls
filt<-het[het$prop_het > 0.20,]
samps<-samps[!samps %in% filt$IID]
##########

## calc pw dist for filtering near identical lines
system("plink --tfile ../data/cowpea --allow-extra-chr --distance square")

### import data and subset to kept lines
d<-read.table("plink.dist")
nm<-read.table("plink.dist.id")
names(d)<-nm$V1
row.names(d)<-nm$V2
d<-d[row.names(d) %in% samps,names(d) %in% samps]

### mtx to pw
xy <- t(combn(colnames(d), 2))
pw<-data.frame(xy, dist=d[xy])

p1<-ggplot(pw, aes(x=dist)) + geom_density() +
	theme_cowplot()
ggsave("dist.pdf", p1, height=4, width=4)

### cluster lines
hc<-hclust(as.dist(d))
p1<-ggdendrogram(hc) +
  geom_hline(aes(yintercept=500), color="red", linetype="dotted")
ggsave("dendro.pdf", p1, height=8, width=8)

### static treecut
grp<-as.data.frame(cutree(hc, h=500))
names(grp)<-"group"
grp$id<-row.names(grp)
write.table(grp, "groups.txt", sep="\t", quote=F, row.names=F)

### check number in each group
tbl<-as.data.frame(table(grp$group))

### plot freq
p1<-ggplot(tbl, aes(x=Freq)) + geom_bar() +
	theme_cowplot()
ggsave("groups.pdf", p1, height=4, width=4)

### sample one id from each group
idlst<-data.frame()
for (i in unique(grp$group)){
	sub<-grp[grp$group==i,]
	set.seed(9)
	sub<-sub[sample(nrow(sub), 1),]
	idlst<-as.data.frame(rbind(idlst, sub))
}

### filter samps
samps<-samps[samps %in% idlst$id]

##########
# write out final datasets
## filter genos
calls<-calls[,names(calls) %in% samps]
df2 <- as.data.frame(unlist(lapply(calls, data.table::tstrsplit, ""),
                              recursive = FALSE))
out<-as.data.frame(cbind(info, df2))
write.table(out, "../data/cowpea.tped", sep=" ", quote=F, row.names=F, col.names=F)

## filter phenos
out2<-out2[out2$V1 %in% samps,]
all.equal(names(calls), out2$V1)
write.table(out2, "../data/cowpea.tfam", sep=" ", quote=F, row.names=F, col.names=F)

## write out sample lst
write.table(samps, "../data/samples.txt", quote=F, row.names=F, col.names=F)
