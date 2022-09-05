#!/usr/bin/env Rscript
# Compare 12-mer abund in TAIR10 vs. Illumina reads
# cjfiscus
# 2022-01-03

setwd("~/Desktop")

library(pacman)
p_load(data.table, ggplot2, cowplot, scico, ggpubr)

df<-fread("~/Desktop/data/ref_vs_illumina_12.txt.gz")
df<-as.data.frame(df)
names(df)<-c("mer", "TAIR10", "Illumina")
df$TAIR10<-log2(df$TAIR10 + 1)
df<-na.omit(df)

## subset for plotting
p1<-ggplot(df, aes(x=TAIR10, y=Illumina)) + stat_bin2d() +
  theme_cowplot() + xlab("12-mer abundance in assembly") +
  ylab("12-mer abundance in Illumina reads") +
  geom_smooth(method="lm", se=FALSE) +
  scale_fill_scico(palette = "lajolla") +
  stat_regline_equation(aes(label=..rr.label..))
ggsave("4_ref_vs_illumina_12mer.pdf", p1, height=5, width=5.5)

# distribution of residuals
## fit model
mod<-lm(Illumina ~ TAIR10, data=df)
res<-as.data.frame(cbind(df$mer, mod$residuals))
names(res)<-c("mer", "residual")

### plot dist of residuals
res$residual<-as.numeric(as.character(res$residual))

p1<-ggplot(res, aes(x=residual)) +
  geom_density() +
  theme_cowplot()
ggsave("5_dist_residual.pdf", p1, height=5, width=5)

## annot rep and non rep
rep<-fread("~/Desktop/data/repetitive_12.txt.gz")
rep$V3<-"repetitive"

nonrep<-fread("~/Desktop/data/non_repetitive_12.txt.gz")
nonrep$V3<-"non-repetitive"

allkmer<-as.data.frame(rbind(rep, nonrep))
allkmer<-allkmer[,c("V1", "V3")]
names(allkmer)<-c("mer", "class")
res1<-merge(res, allkmer, by="mer")

### plot distribution of residuals per class
#### cleanup
rm(allkmer);rm(rep);rm(nonrep)
res1$residual<-as.numeric(as.character(res1$residual))

cnt<-as.data.frame(table(res1$mer))
cnt<-cnt[cnt$Freq == 1,]
s<-res1[res1$mer %in% cnt$Var1,]

## raw residuals
p1<-ggplot(s, aes(x=class, y=residual)) +
  geom_violin() +
  theme_cowplot() +
  stat_summary(fun="median", geom="point") +
  geom_hline(aes(yintercept=0), linetype="dotted")
  
ggsave("6_dist_residual_class.pdf", p1, height=5, width=5)

## abs residuals
p1<-ggplot(s, aes(x=class, y=abs(residual))) +
  geom_violin() +
  theme_cowplot() +
  stat_summary(fun="median", geom="point") +
  geom_hline(aes(yintercept=0), linetype="dotted")
ggsave("6_dist_residual_class_abs.pdf", p1, height=5, width=5)
