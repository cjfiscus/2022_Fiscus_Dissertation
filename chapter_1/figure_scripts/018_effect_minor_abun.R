#!/usr/bin/env Rscript
# determine effect of minor alleles
# cjfiscus
# 2021-09-08

library(pacman)
p_load(data.table, ggplot2, cowplot, dplyr, ggsci)

# import list of snps to consider
df<-fread("~/Desktop/working/arabidopsis_kmer_final_time/metagwas_final/metagwas_top1326.txt")

# import lst of effects
retro<-fread("~/Desktop/working/arabidopsis_kmer_final_time/metagwas_final/retrotrans_beta_mean.txt")
retro$class<-"Retrotransposon"
dna<-fread("~/Desktop/working/arabidopsis_kmer_final_time/metagwas_final/dnatrans_beta_mean.txt")
dna$class<-"DNA transposon"
sr<-fread("~/Desktop/working/arabidopsis_kmer_final_time/metagwas_final/sr_beta_mean.txt")
sr$class<-"Simple Repeat"

## combine into one df
effects<-as.data.frame(rbind(retro, dna, sr))

## cleanup
rm(retro);rm(dna);rm(sr)

## merge data
m<-merge(df, effects, by=c("snp_id", "class"))
m$effect<-ifelse(m$beta_mean > 0, "increase", "decrease")

## plot for variants considered individually
m$effect<-factor(m$effect, levels=c("increase", "decrease"))
m$class1<-factor(m$class, levels=c("Retrotransposon", "DNA transposon", "Simple Repeat"))
p1<-ggplot(m, aes(x=effect)) + geom_bar() +
  facet_grid(.~class1) +
  theme_cowplot() + xlab("minor allele effect on copy number")
ggsave("~/Desktop/18_effect_by_class.pdf", p1, height=2, width=6)

## subset to snps found in two or more classes
tbl<-as.data.frame(table(m$snp_id))
tbl<-tbl[tbl$Freq > 1,]

m1<-m[m$snp_id %in% tbl$Var1,]

## consider pairwise
lst<-data.frame()
for (i in unique(m1$snp_id)){
  sub<-m1[m1$snp_id == i,]
  sub<-sub[,c("snp_id", "class", "effect")]
  
  sub1<-sub
  names(sub1)<-c("snp_id", "class1", "effect1")
  
  sub2<-sub
  names(sub2)<-c("snp_id", "class2", "effect2")
  
  combs<-as.data.frame(t(combn(sub$class, 2)))
  names(combs)<-c("class1", "class2")
  
  sub3<-merge(combs, sub1, by=c("class1"))
  sub3<-merge(sub3, sub2, by=c("class2"))
  
  sub3<-sub3[,c("snp_id.x", "class1", "class2", "effect1", "effect2")]
  lst<-as.data.frame(rbind(lst, sub3))
  
}

## process lst for plotting
names(lst)[1]<-"snp_id"
lst$effect<-ifelse(lst$effect1=="increase" & lst$effect2 == "increase" | 
                     lst$effect1 == "decrease" & lst$effect2 == "decrease", 
                   as.character(lst$effect1), paste0(lst$effect1, "/",lst$effect2))

lst$class <-paste0(lst$class1, " & ",lst$class2 )

lst$effect_dir<-ifelse(lst$effect=="increase" | lst$effect =="decrease", "nonantagonistic", "antagonistic")
lst$effect_dir<-factor(lst$effect_dir, levels=c("nonantagonistic", "antagonistic"))

lst2<-as.data.frame(table(lst$class, lst$effect_dir, lst$effect))
ag<-as.data.frame(aggregate(Freq ~ Var1, lst2, FUN="sum"))
names(ag)[2]<-"s"

lst2<-merge(lst2, ag, by="Var1")
lst2$prop<-as.numeric(lst2$Freq/lst2$s)
lst2<-lst2[lst2$prop > 0,]

## make plots
lst2$Var3<-factor(lst2$Var3, levels=c("increase", "decrease", "increase/decrease",
                                      "decrease/increase"))
names(lst2)[3]<-"direction"
p1<-ggplot(lst2, aes(x=Var2, y=prop, fill=direction)) + geom_bar(stat="identity") + 
  facet_grid(.~Var1) + theme_cowplot() + theme(legend.position="bottom") +
  ylab("Proportion of sites") + xlab("effect") + ylim(0, 0.6) +
  scale_fill_npg()
ggsave("~/Desktop/18_effect_by_class_pw.pdf", p1, height=3, width=10)
