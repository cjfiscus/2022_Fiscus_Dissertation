#!/usr/bin/env Rscript
# Figure 3
# cjfiscus
# 2022-04-22

library(pacman)
p_load(ggplot2, cowplot, ggsci, tidyr, stringr, readxl, 
       patchwork, dplyr, plyr, ggpubr, reshape2, ggrepel)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

## 3A phenotypes into categories with tests
# import data
df<-read.delim("../DATA/master_phenotype_table.txt")

df1<-read.delim("../DATA/pheno_by_group_wilcox.txt")
df1$sig<-ifelse(df1$wilcox_pvalue < 0.05/nrow(df1), "significant", "not significant")
df1<-df1[,c("ID", "sig")]
df1$q<-"quantitative"

df2<-read.delim("../DATA/pheno_by_group_fisher.txt")
df2$sig<-ifelse(df2$fisher_pvalue < 0.05/nrow(df2), "significant", "not significant")
df2<-df2[,c("ID", "sig")]
df2$q<-"qualitative"
df1<-as.data.frame(rbind(df1, df2))
df<-merge(df, df1, by="ID")

tbl<-as.data.frame(table(df$CATEGORY_2, df$sig, df$q))

tbl1<-as.data.frame(table(df$CATEGORY_2))
tbl1<-tbl1[order(tbl1$Freq),]
tbl$Var1<-factor(tbl$Var1, levels=tbl1$Var1)

#tbl$Var2<-factor(tbl$Var2, levels=c("significant", "not significant"))
tbl$class<-paste0(tbl$Var2, "_", tbl$Var3)

pal1<-c("dimgrey", "darkgoldenrod1")
names(pal1)<-c("not significant", "significant")

# plot for quantitative traits
tbl$Var3<-factor(tbl$Var3, levels=c("quantitative", "qualitative"))
tbl$Var2<-factor(tbl$Var2, levels=c("significant", "not significant"))
p1<-ggplot(tbl, aes(x=Freq, y=Var1, fill=Var2)) + geom_bar(stat="identity") +
  theme_cowplot() + ylab("Category") + xlab("Frequency") +
  theme(legend.position="none", strip.background = element_blank()) +
  scale_fill_manual(values=pal1) +
  facet_grid(.~Var3, scales="free_x")

##########
# quant phenotype boxplots
# import data
df<-read.delim("../data/cluster_assignments.txt")

df1<-read_excel("../data/IITA-Cowpea collection.xls")
df1<-df1[,c("Accession name", "Latitude", "Longitude", "Country of origin", "Biological status of accession")]
names(df1)[1]<-"sample"

## merge datasets
df<-merge(df, df1, by="sample")
k3<-df[df$K==2,]
k3<-as.data.frame(unique(k3[,c("sample", "assignment")]))

# import phenotypic dataset
## full list of phenos to consider
phenos<-read.delim("../data/master_phenotype_table.txt")

### combined phenotypes
df<-read.table("../data/cowpea_gwas_combined.fam")
nm<-read.table("../data/combined_phenotypes.txt", header=F, sep = "\t")
names(df)[6:ncol(df)]<-nm$V2
df<-df[,c(2, 6:ncol(df))]
names(df)[1]<-"sample"
w<-melt(df, id.vars=c("sample"))

### uncombined phenotypes
df<-read.table("../data/cowpea_gwas.fam")
nm<-read.table("../data/phenotypes.txt", header=F, sep="\t")
names(df)[6:ncol(df)]<-nm$V2
df<-df[,c(2, 6:ncol(df))]
names(df)[1]<-"sample"
df[,2:16]<-(exp(1) ^ df[,2:16])-0.01
w1<-melt(df, id.vars=c("sample"))

### merge
w<-as.data.frame(rbind(w, w1))

### subset to phenos we studied
w<-w[w$variable %in% phenos$PHENOTYPE,]

## wide to long
w<-merge(w, k3, by="sample")

## plots for sig quant pheno
names(w)[2]<-"PHENOTYPE"
w<-merge(w, phenos[,c("ID", "PHENOTYPE")], by="PHENOTYPE")
w$assignment<-factor(w$assignment, levels=c("Cluster1", "Admixed", "Cluster2"))

### subset to sig
df1<-read.delim("../DATA/pheno_by_group_wilcox.txt")
a<-0.05/nrow(df1)
df1<-df1[df1$wilcox_pvalue < a,]

for (i in unique(df1$ID)){
  sub<-w[w$ID==i,]
  
  ## plotting
  pal2<-c("#D55E00", "#009E73", "gray69")
  names(pal2)<-c("Cluster1", "Cluster2", "Admixed")
  cnt<-as.data.frame(table(na.omit(sub, na.rm=T)$assignment))
  #cnt$Freq<-paste0("N = ",cnt$Freq)
  names(cnt)<-c("assignment", "N")
  sub<-merge(sub, cnt, by="assignment")
  comps<-list(c(as.character(cnt[1,2]), as.character(cnt[3,2])))
  sub$N<-factor(sub$N, levels=cnt$N)
  
  p<-ggplot(sub, aes(x=N, y=value)) + 
    geom_boxplot(aes(color=assignment)) +
    theme_cowplot() +
    stat_compare_means(comparisons=comps, label = "p.signif") +
    xlab("group") +
    ylab(as.character(unique(sub$PHENOTYPE))) +
    scale_color_manual(values=pal2) + 
    theme(legend.position="none",
          axis.title.x= element_blank(),
          axis.ticks.x=element_blank())
  assign(paste0("b", i), p)
}
##########
# plot frequency per group
df<-read.delim("../DATA/pheno_by_group_fisher.txt")
df$sig<-ifelse(df$fisher_pvalue < 0.05/nrow(df), "significant", "not significant")
df$c1c2<-abs(df$freq1_Cluster1-df$freq1_Cluster2)
df$c1a<-abs(df$freq1_Cluster1-df$freq1_Admixed)
df$c2a<-abs(df$freq1_Cluster2-df$freq1_Admixed)

pal1<-c("dimgrey", "darkgoldenrod1")
names(pal1)<-c("not significant", "significant")

## cluster 1 vs cluster 2
p2<-ggplot(df, aes(x=freq1_Cluster1, y=freq1_Cluster2)) +
  geom_abline(intercept = 0, alpha=0.5) +
  geom_point(aes(color=sig), alpha=0.7) +
  geom_text_repel(data=df[df$c1c2 >=0.10,], 
                  aes(x=freq1_Cluster1, y=freq1_Cluster2, label=ID),
                  box.padding = 0.5, max.overlaps = Inf) +
  theme_cowplot() +
  theme(legend.position="none") +
  xlab("Cluster 1 frequency") +
  ylab("Cluster 2 frequency") +
  scale_color_manual(values=pal1)

## cluster 1 vs admixed
p3<-ggplot(df, aes(x=freq1_Cluster1, y=freq1_Admixed)) +
  geom_abline(intercept = 0, alpha=0.5) +
  geom_point(aes(color=sig), alpha=0.7, color="dimgrey") +
  geom_text_repel(data=df[df$c1a >=0.10,], 
                  aes(x=freq1_Cluster1, y=freq1_Admixed, label=ID),
                  box.padding = 0.5, max.overlaps = Inf) +
  theme_cowplot() +
  theme(legend.position="none") +
  xlab("Cluster 1 frequency") +
  ylab("Admixed frequency")

## cluster 2 vs admixed
p4<-ggplot(df, aes(x=freq1_Cluster2, y=freq1_Admixed)) +
  geom_abline(intercept = 0, alpha=0.5) +
  geom_point(aes(color=sig), alpha=0.7, color="dimgrey") +
  geom_text_repel(data=df[df$c2a >=0.10,], 
                  aes(x=freq1_Cluster2, y=freq1_Admixed, label=ID),
                  box.padding = 0.5, max.overlaps = Inf) +
  theme_cowplot() +
  theme(legend.position="none") +
  xlab("Cluster 2 frequency") +
  ylab("Admixed frequency")

##########
# make main figure
layout<-"
AAAAAAAA
BCDEFGHI
JJKKLL##
"
row1<-p1 + b1 + b2 + plot_layout(nrow=1, widths=c(4,1,1))

row1a<-p1
ggsave("~/Desktop/row1a.pdf", p1, height=4, width=8)

row2<-b4 + 
  b5 + b6 + b7 + b8 + 
  b9 + plot_layout(nrow=1)
row3<-p2 + p3 + p4

fig2<-row2/row3
ggsave("~/Desktop/row2_3.pdf", fig2, height=8, width=12)

##########
# make sup figure
lst<-c("4C","5C", "9C", "28C", "29C", "30C", "31C", "32C", "40C")

for (i in lst){
  sub<-w[w$ID==i,]
  pal2<-c("#D55E00", "#009E73", "gray69")
  names(pal2)<-c("Cluster1", "Cluster2", "Admixed")
  cnt<-as.data.frame(table(na.omit(sub, na.rm=T)$assignment))
  names(cnt)<-c("assignment", "N")
  sub<-merge(sub, cnt, by="assignment")
  comps<-list(c(as.character(cnt[1,2]), as.character(cnt[3,2])))
  sub$N<-factor(sub$N, levels=cnt$N)
  
  tbl<-as.data.frame(table(sub$assignment, sub$value))
  agg<-as.data.frame(aggregate(Freq ~ Var1, data=tbl, FUN=sum))
  names(agg)[2]<-"total"
  tbl<-merge(tbl, agg, by="Var1")
  tbl$prop<-tbl$Freq/tbl$total
  
  p<-ggplot(tbl, aes(x=Var2, y=prop)) + 
    geom_bar(aes(fill=Var1), stat="identity", position="dodge") +
    theme_cowplot() +
    stat_compare_means(comparisons=comps, label = "p.signif") +
    xlab("value") +
    ylab("Proportion samples per group") +
    scale_fill_manual(values=pal2) + 
    theme(legend.position="none") +
    ggtitle(as.character(unique(sub$PHENOTYPE)))
  assign(paste0("b", i), p)
}

row<-b28C + b29C + b30C + b31C + b32C + b40C + b4C + b5C + b9C +
  plot_layout(nrow=3, ncol=4)
ggsave("~/Desktop/supp.pdf", row, height=9, width=12)

leg<-get_legend(b28C + theme(legend.position="right"))
ggsave("~/Desktop/fig2leg.pdf", leg, height=3, width=3)

