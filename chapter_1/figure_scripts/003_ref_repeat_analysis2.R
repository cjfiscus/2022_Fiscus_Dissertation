3 #!/usr/bin/env Rscriptf
# ref repeat analysis
# cjfiscus
# 2022-01-05

setwd("~/Desktop")

library(pacman)
p_load(ggplot2, cowplot, ggpubr)

# 12-mer based seq abund ests
## norm seq abund ests
df<-read.table("TAIR10_seq_abun_TE.txt")
names(df)<-c("id", "abun", "N")

## parse TE families
df$id2<-sapply(strsplit(df$id, split="\\|"), "[", 5)

df1<-read.table("TAIR10_seq_abun_buscos.txt")
names(df1)<-c("id", "abun", "N")

### calc med of buscos and use for norm
m<-median(df1$abun)
df$norm<-df$abun/m

### for LTR retro use internal seqs
dup<-df
dup$id2<-paste0(df$id2, "_I")
df<-as.data.frame(rbind(df, dup))
df<-df[,c("id2", "norm", "N")]
names(df)[1]<-"ID"

### account for proportion of seq in genome
rpt<-read.delim("./data/A_thaliana_Repeat_Database.txt")
rpt<-rpt[,c("ID", "LENGTH", "CLASS", "SUBCLASS")]
df<-merge(df, rpt, by="ID")
df$scalar<-df$N/df$LENGTH
df$norm<-df$norm*df$scalar

## tally abund by family
df<-as.data.frame(aggregate(norm ~ ID, data=df, FUN=sum))
names(df)<-c("id", "kmer")

##########
## annot based seq abun ests
annot<-read.delim("TAIR10_Transposable_Elements.txt")
annot$len<-annot$Transposon_max_End-annot$Transposon_min_Start
annot<-annot[,c("Transposon_Name", "Transposon_Family", "len")]
#annot<-annot[,c("Transposon_Name", "Transposon_Family")]
names(annot)[2]<-"ID"

rpt<-read.delim("./data/A_thaliana_Repeat_Database.txt")
rpt<-rpt[,c("ID", "LENGTH", "CLASS", "SUBCLASS")]

### use internal seqs for LTR retrotrans
ltr<-annot[!annot$ID %in% rpt$ID,]
ltr$ID<-paste0(ltr$ID, "_I")
annot<-as.data.frame(rbind(annot, ltr))

annot<-merge(annot, rpt, by="ID")
annot$cn<-annot$len/annot$LENGTH

#### count 
df2<-as.data.frame(aggregate(cn ~ ID, data=annot, FUN=sum))
names(df2)<-c("id", "annot")

# compare
m<-merge(df, df2, by="id")

p1<-ggplot(m, aes(x=annot, y=kmer)) + 
  geom_point() +
  theme_cowplot() +
  geom_smooth(method="lm") +
  stat_regline_equation(aes(label=..rr.label..)) +
  xlab("annotation based estimate") +
  ylab("12-mer based estimate")
ggsave("3_annot_12mer_ref.pdf", p1, height=4, width=4)
##########

# compare these estimates to 12-mer based estimates from Illumina reads
df<-read.delim("~/Desktop/data/A_thal_kmer_abund_norm.txt", check.names=F)
df<-df[,c("Feature", "6909")]
names(df)<-c("id", "illumina")

m1<-merge(m, df, by="id")

p1<-ggplot(m1, aes(x=kmer, y=illumina)) +
  geom_point() +
  theme_cowplot() +
  geom_smooth(method="lm") +
  stat_regline_equation(aes(label=..rr.label..)) +
  xlab("annotation based estimate (reference)") +
  ylab("12-mer based estimate (Illumina)")
ggsave("8_ref_illumina_12mer_est.pdf", p1, height=4, width=4)

p1<-ggplot(m1, aes(x=log(kmer), y=log(illumina))) +
  geom_point() +
  theme_cowplot() +
  geom_smooth(method="lm") +
  stat_regline_equation(aes(label=..rr.label..)) +
  xlab("annotation based estimate (reference)") +
  ylab("12-mer based estimate (Illumina)")
ggsave("8_ref_illumina_12mer_est_loglog.pdf", p1, height=4, width=4)
