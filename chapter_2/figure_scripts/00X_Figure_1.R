#!/usr/bin/env Rscript
# Figure 1
# cjfiscus
# 2022-07-22

library(pacman)
p_load(ggplot2, cowplot, ggridges, scales, patchwork, ggforce)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/SCRIPTS/")

## first panel scans of variation compared to gene and repeat content
df<-read.table("../data/cbpco_feats_100kb.txt", header=T)
df1<-read.table("../data/cbpcr_feats_100kb.txt", header=T)
df<-as.data.frame(rbind(df, df1))
df$chr<-factor(df$chr, levels=paste0("SCF_", 1:16))
df$Var2<-factor(df$Var2, levels=c("gene", "repeats", "insertion", "deletion",
                                  "translocation", "inversion"))
df$subg<-ifelse(df$chr %in% paste0("SCF_", 1:8), "CbpCo", "CbpCr")
df$chr2<-sub("SCF_", "", df$chr)
df$chr2<-factor(df$chr2, levels=1:16)

pal<-c("#7570b3", "#1b9e77")
names(pal)<-c("CbpCo", "CbpCr")

p1<-ggplot(df, aes(x=Var1, y=Freq, fill=subg)) + geom_bar(stat="identity") + 
  theme_cowplot() +
  scale_y_continuous(breaks=pretty_breaks()) +
  facet_grid(Var2~chr2, scales="free") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.title.x=element_blank(), 
        legend.position="none") +
  scale_fill_manual(values=pal)
ggsave("~/Desktop/1A.pdf", p1, height=4, width=12)
##########

## feature counts
tbl<-as.data.frame(aggregate(Freq ~ subg + Var2, data=df, FUN = sum))
tbl<-tbl[!tbl$Var2=="gene",]
tbl<-tbl[!tbl$Var2=="repeats",]

p2<-ggplot(tbl, aes(x=Var2, y=Freq, fill=subg)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=pal) + theme_cowplot() + 
  theme(legend.position="none", axis.title.x = element_blank()) + xlab("SV class")
##########
## feature lengths
df<-read.delim("../data/Co_CbpCo_sv.txt.gz")
df$spp<-"CbpCo"
df1<-read.delim("../data/Cr_CbpCr_sv.txt.gz")
df1$spp<-"CbpCr"
df<-as.data.frame(rbind(df, df1))
df$type<-factor(df$type, levels=c("insertion", "deletion", "translocation",
                                  "inversion"))
df$spp<-factor(df$spp, levels=c("CbpCo", "CbpCr"))
df$len<-abs(df$length)+1
df$len<-as.numeric(df$len)

p3<-ggplot(df, aes(x=len, y=spp, fill=spp)) + 
  geom_density_ridges(alpha=0.5) + facet_grid(.~type) + 
  theme_cowplot() + scale_fill_manual(values=pal) +
  theme(legend.position="none", strip.background = element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(trans="log", breaks=log_breaks()) +
  xlab("size (bp)")
##########
## SNP frequency plot
df<-read.delim("../data/Co_CbpCo_snp.txt.gz")
df$spp<-"CbpCo"
df1<-read.delim("../data/Cr_CbpCr_snp.txt.gz")
df1$spp<-"CbpCr"
df<-as.data.frame(rbind(df, df1))
df$proportion<-df$Freq/df$total

agg<-as.data.frame(aggregate(proportion ~ Var1, data=df, FUN=mean))
agg<-agg[order(agg$proportion, decreasing = T),]
df$Var1<-factor(df$Var1, levels=agg$Var1)

p4<-ggplot(df, aes(x=Var1, y=proportion, fill=spp)) + geom_bar(stat="identity", 
                                                           position="dodge") +
  xlab("substitution") + theme_cowplot() + theme(legend.position="top") +
  scale_fill_manual(values=pal) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())
##########
# pg plots
### overview
df<-read.delim("../data/pg_overall.txt", header=T)
df$cat<-factor(df$cat, levels=c("no CNV", "CNV"))

p5<-ggplot(df, aes(x=cat, y=freq)) + geom_bar(stat="identity") +
  theme_cowplot() + theme(axis.title.x = element_blank())

### increase
df<-read.delim("../data/pg_shared.txt")
df$spp<-factor(df$spp, levels=c("CbpCo", "CbpCr", "both"))

df1<-df[df$cat=="increase",]

pal2<-c("#7570b3", "#1b9e77", "#488795")
names(pal2)<-c("CbpCo", "CbpCr", "both")

p6<-ggplot(df1, aes(x=spp, y=Freq, fill=spp)) + geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_manual(values=pal2) + theme(legend.position="none", 
                                         axis.title.x = element_blank()) +
  ggtitle("increase")

##### PLOT FOR PRESENTATION
df$cat<-factor(df$cat, levels=c("no change", "increase", "decrease", "lost"))

p00<-ggplot(df, aes(x=cat, y=Freq, fill=spp)) + 
  geom_bar(stat="identity", position = position_dodge2(width = 1, preserve = "single")) +
  theme_cowplot() +
  scale_fill_manual(values=pal2) + theme(legend.position="none", 
                                         axis.title.x = element_blank())
ggsave("~/Desktop/genes.pdf", p00, height=4, width=5)
leg<-get_legend(p00 + theme(legend.position="top"))
ggsave("~/Desktop/leg2.pdf", leg, height=2, width=5)
#####

### decrease
df1<-df[df$cat=="decrease",]
p7<-ggplot(df1, aes(x=spp, y=Freq, fill=spp)) + geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_manual(values=pal2) + theme(legend.position="none", 
                                         axis.title.x = element_blank()) +
  ggtitle("decrease")

### lost
df1<-df[df$cat=="lost",]
p8<-ggplot(df1, aes(x=spp, y=Freq, fill=spp)) + geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_manual(values=pal2) + theme(legend.position="none", 
                                         axis.title.x = element_blank()) +
  ggtitle("lost")
##########
# plot GO terms
## increase
df<-read.delim("../data/cn_increase_GO.txt")
df<-df[1:10,]
df$Term<-factor(df$Term, levels=df$Term)
df$enrich<-as.numeric(df$enrich)

p9<-ggplot(df, aes(x=Term, y=enrich, group=Term)) + geom_point() +
  theme_cowplot() + theme(axis.text.x=element_text(angle = 45, hjust=1), 
                          axis.title.x=element_blank()) +
  ylab("Enrichment") + 
  scale_y_continuous(breaks=pretty_breaks(), limits = c(0, ceiling(max(df$enrich))))

p9a<-ggplot(df, aes(x=Term, y=enrich, group=Term)) + geom_point() +
  theme_cowplot() + theme(axis.text.x=element_blank(), 
                          axis.title.x=element_blank()) +
  ylab("Enrichment") + 
  scale_y_continuous(breaks=pretty_breaks(), limits = c(0, ceiling(max(df$enrich))))

## lost
df<-read.delim("../data/cn_lostall_GO.txt")
df<-df[1:10,]
df$Term<-factor(df$Term, levels=df$Term)
df$enrich<-as.numeric(df$enrich)

p10<-ggplot(df, aes(x=Term, y=enrich, group=Term)) + geom_point() +
  theme_cowplot() + theme(axis.text.x=element_text(angle = 45, hjust=1), 
                          axis.title.x=element_blank()) +
  ylab("Enrichment") + 
  scale_y_continuous(breaks=pretty_breaks(), limits = c(0, ceiling(max(df$enrich))))

p10a<-ggplot(df, aes(x=Term, y=enrich, group=Term)) + geom_point() +
  theme_cowplot() + theme(axis.text.x=element_blank(), 
                          axis.title.x=element_blank()) +
  ylab("Enrichment") + 
  scale_y_continuous(breaks=pretty_breaks(), limits = c(0, ceiling(max(df$enrich))))

##########
# plot lineage specific change
df<-read.delim("../data/pg_single_lineage.txt")
df<-df[!df$cat=="unique",]
df$cat<-factor(df$cat, levels=c("no CNV", "increase", "decrease", "lost",
                                "translocated"))
df1<-df[df$donor=="Co",]

p11<-ggplot(df1, aes(x=cat, y=Freq, fill=spp)) + geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_manual(values=pal) +
  theme(legend.position="none", axis.title.x = element_blank())
df1<-df[df$donor=="Cr",]
p12<-ggplot(df1, aes(x=cat, y=Freq, fill=spp)) + geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_manual(values=pal) +
  theme(legend.position="none", axis.title.x = element_blank())
###########
# plot repeat plots
### copy number
df<-read.table("../data/repeat_summary.txt")
names(df)<-c("spp", "class", "family", "cn")
keep<-c("CbpCo", "CbpCr", "Co39", "Cr145")
df<-df[df$spp %in% keep,]

pal3<-c("#ACA9BB", "#7570b3", "#1b9e77", "#C8FCEA")
names(pal3)<-c("Co39", "CbpCo", "CbpCr", "Cr145")
df$spp<-factor(df$spp, levels=c("Co39", "CbpCo", "CbpCr", "Cr145"))
df$class<-factor(df$class, levels=c("LTR", "nonLTR", "TIR", "nonTIR"))

p13<-ggplot(df, aes(x=spp, y=cn, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~class, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position="none") +
  ylab("copy number")

### by family (supplement)
df1<-df[df$class=="LTR",]

p13a<-ggplot(df1, aes(x=spp, y=cn, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position="none") +
  ylab("copy number")

df1<-df[df$class=="nonLTR",]

p13b<-ggplot(df1, aes(x=spp, y=cn, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position="none") +
  ylab("copy number")

df1<-df[df$class=="TIR",]

p13c<-ggplot(df1, aes(x=spp, y=cn, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position="none") +
  ylab("copy number")

df1<-df[df$class=="nonTIR",]

p13d<-ggplot(df1, aes(x=spp, y=cn, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position="none") +
  ylab("copy number")

plot<-ggplot(df, aes(x=family, y=cn, fill=spp)) + 
  geom_bar(stat="identity", position="dodge") +
  facet_grid(.~class, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(), legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("copy number") + facet_grid(.~class, scales="free", space="free")
  
ggsave("~/Desktop/S2.pdf", plot, height=4, width=12)

### sequence
df<-read.table("../data/repeat_summary2.txt")

names(df)<-c("spp", "class", "family", "bp")
keep<-c("CbpCo", "CbpCr", "Co39", "Cr145")
df<-df[df$spp %in% keep,]

pal3<-c("#ACA9BB", "#7570b3", "#1b9e77", "#C8FCEA")
names(pal3)<-c("Co39", "CbpCo", "CbpCr", "Cr145")
df$spp<-factor(df$spp, levels=c("Co39", "CbpCo", "CbpCr", "Cr145"))
df$class<-factor(df$class, levels=c("LTR", "nonLTR", "TIR", "nonTIR"))

p14<-ggplot(df, aes(x=spp, y=bp, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~class, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position="none") +
  ylab("sequence (Mbp)") +
  scale_y_continuous(labels=function(x)x/1000000)

### by family (supplement)
df1<-df[df$class=="LTR",]

p14a<-ggplot(df1, aes(x=spp, y=bp, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="none") +
  ylab("sequence (Mbp)") +
  scale_y_continuous(labels=function(x)x/1000000)

df1<-df[df$class=="nonLTR",]

p14b<-ggplot(df1, aes(x=spp, y=bp, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="none") +
  ylab("sequence (Mbp)") +
  scale_y_continuous(labels=function(x)x/1000000)

df1<-df[df$class=="TIR",]

p14c<-ggplot(df1, aes(x=spp, y=bp, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="none") +
  ylab("sequence (Mbp)") +
  scale_y_continuous(labels=function(x)x/1000000)

df1<-df[df$class=="nonTIR",]

p14d<-ggplot(df1, aes(x=spp, y=bp, fill=spp)) + 
  geom_bar(stat="identity") +
  facet_grid(.~family, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        legend.position="none") +
  ylab("sequence (Mbp)") +
  scale_y_continuous(labels=function(x)x/1000000)
leg<-get_legend(p13 + theme(legend.position="right"))

plot<-ggplot(df, aes(x=family, y=bp, fill=spp)) + 
  geom_bar(stat="identity", position="dodge") +
  facet_grid(.~class, scales = "free") + theme_cowplot() +
  scale_fill_manual(values=pal3) +
  theme(strip.background = element_blank(), legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("sequence (Mbp)") + facet_grid(.~class, scales="free", space="free") +
  scale_y_continuous(labels=function(x)x/1000000)

ggsave("~/Desktop/S3.pdf", plot, height=4, width=12)
ggsave("~/Desktop/leg.pdf", leg, height=2, width=2)

###########
## write out components of figure
ggsave("~/Desktop/A.pdf", p1, height=3, width=10)
ggsave("~/Desktop/B.pdf", p2, height=3, width=3)
ggsave("~/Desktop/C.pdf", p3, height=3, width=4)
ggsave("~/Desktop/D.pdf", p4, height=3, width=4)
ggsave("~/Desktop/E.pdf", p6, height=3, width=1.5)
ggsave("~/Desktop/F.pdf", p9, height=4, width=3)
ggsave("~/Desktop/F1.pdf", p9a, height=1.5, width=3)
ggsave("~/Desktop/G.pdf", p7, height=3, width=1.5)
ggsave("~/Desktop/H.pdf", p8, height=3, width=1.5)
ggsave("~/Desktop/I.pdf", p10, height=4, width=3)
ggsave("~/Desktop/I1.pdf", p10a, height=1.5, width=3)
ggsave("~/Desktop/J.pdf", p11, height=3, width=3)
ggsave("~/Desktop/K.pdf", p12, height=3, width=3)
ggsave("~/Desktop/L.pdf", p13, height=3, width=5)
ggsave("~/Desktop/M.pdf", p14, height=3, width=5)
ggsave("~/Desktop/leg.pdf", leg, height=2, width=3)

### repeat subplots 
ggsave("~/Desktop/13a.pdf", p13a, height=3, width=8)
ggsave("~/Desktop/13b.pdf", p13b, height=3, width=8)
ggsave("~/Desktop/13c.pdf", p13c, height=3, width=8)
ggsave("~/Desktop/13d.pdf", p13d, height=3, width=8)

ggsave("~/Desktop/14a.pdf", p14a, height=3, width=8)
ggsave("~/Desktop/14b.pdf", p14b, height=3, width=8)
ggsave("~/Desktop/14c.pdf", p14c, height=3, width=8)
ggsave("~/Desktop/14d.pdf", p14d, height=3, width=8)
