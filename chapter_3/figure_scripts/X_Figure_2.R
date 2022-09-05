#!/usr/bin/env Rscript
# Figure 1
# cjfiscus
# 2022-04-08

library(pacman)
p_load(ggplot2, cowplot, ggsci, tidyr, stringr, readxl, 
       patchwork, dplyr, plyr, ggpubr)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

## 1A delta K
# import data
df<-read.table("../DATA/logli.txt")
names(df)<-c("K", "logli")

## calc mean and sd of logli
d<-as.data.frame(aggregate(logli ~ K, data=df, FUN=mean))
names(d)<-c("K", "mlogli")
e<-as.data.frame(aggregate(logli ~ K, data=df, FUN=sd))
names(e)<-c("K", "slogli")
d<-merge(d, e, by="K")
d$delta<-d$mlogli-lag(d$mlogli)
d$delta2<-abs(lead(d$delta)-d$delta)
d$deltak<-d$delta2/d$slogli

dK<-ggplot(d, aes(x=K, y=deltak)) + geom_line() + geom_point() +
  theme_cowplot() +
  scale_x_continuous(breaks=seq(2,9,1), limits=c(2,9)) +
  ylab(expression(Delta*"K")) +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank())
##########
##########
## admixture barplots
# pal for this section
pal<-c("#D55E00", "#009E73", "#56B4E9", "#CC79A7")
names(pal)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4")

### k = 2
df<-read.delim("../data/cluster_assignments.txt")
k2<-df[df$K==2,]

## produce sort order
ord<-k2[,c("sample", "value", "variable")]
ord<-spread(ord, variable, value)
maxval<-apply(ord[,2:ncol(ord)], 1, max)
matchval<-vector(length=nrow(ord))
for(j in 1:nrow(ord)) matchval[j] <- match(maxval[j], ord[j, ])
ord$maxval<-maxval
ord$matchval<-matchval
ord<-ord[with(ord, order(matchval, -maxval)),]

k2$sample<-factor(k2$sample, levels=ord$sample)

## add regional data
samps<-read_excel("../data/IITA-Cowpea collection.xls")
samps<-samps[samps$`Accession name` %in% k2$sample,]
samps<-samps[,c("Accession name", "Biological status of accession", 
                "Country of origin", "Collecting/acquisition source", "Latitude",
                "Longitude")]
unsd<-read.csv("../data/UNSD_M49_2022.csv")
unsd<-unsd[,c("Country.or.Area", "Region.Name", "Sub.region.Name", "Intermediate.Region.Name")]
names(unsd)[1]<-"Country of origin"

## reval countries so datasets match
samps$`Country of origin`<-revalue(samps$`Country of origin`, 
                                   c("Cote d'Ivoire"="Côte d’Ivoire",
                                     "Iran"="Iran (Islamic Republic of)",
                                     "Republic of the Congo"="Democratic Republic of the Congo",
                                     "Tanzania"="United Republic of Tanzania",
                                     "United Kingdom"="United Kingdom of Great Britain and Northern Ireland",
                                     "Zaire"="Democratic Republic of the Congo",
                                     "Swaziland"="Eswatini",
                                     "Union of Soviet Socialist Republics"="Ukraine",
                                     "Venezuela"="Venezuela (Bolivarian Republic of)"))

samps<-as.data.frame(unique(merge(samps, unsd, by="Country of origin", all.x=T)))
names(samps)[2]<-"sample"
k2<-merge(k2, samps, by="sample", all.x=T)

sub1<-k2
sub1$Intermediate.Region.Name<-ifelse(sub1$Sub.region.Name=="Northern Africa", "Northern", sub1$Intermediate.Region.Name)
sub1$Intermediate.Region.Name<-revalue(sub1$Intermediate.Region.Name, c("Eastern Africa"="Eastern", "Middle Africa"="Middle", 
                                                                        "Southern Africa"="Southern", "Western Africa"="Western"))
sub1$Intermediate.Region.Name<-factor(sub1$Intermediate.Region.Name, levels=c("Northern","Western", "Middle", "Eastern", "Southern"))
sub1<-sub1[complete.cases(sub1$`Country of origin`),]

k2p<-ggplot(sub1, aes(x=sample, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  theme_cowplot() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(), 
        strip.background = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        axis.line=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  ylab("Proportion of ancestry") +
  xlab("accession") +
  scale_fill_manual(values=pal) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(.~Region.Name+Intermediate.Region.Name, scales="free_x", space="free_x")

ggsave("~/Desktop/k2p.pdf", k2p, height=3, width=10)
#########
### k = 3
# import data
k3<-df[df$K==3,]

## produce sort order
ord<-k3[,c("sample", "value", "variable")]
ord<-spread(ord, variable, value)
maxval<-apply(ord[,2:ncol(ord)], 1, max)
matchval<-vector(length=nrow(ord))
for(j in 1:nrow(ord)) matchval[j] <- match(maxval[j], ord[j, ])
ord$maxval<-maxval
ord$matchval<-matchval
ord<-ord[with(ord, order(matchval, -maxval)),]

k3$sample<-factor(k3$sample, levels=ord$sample)
k3<-merge(k3, samps, by="sample", all.x=T)

sub1<-k3
sub1$Intermediate.Region.Name<-ifelse(sub1$Sub.region.Name=="Northern Africa", "Northern", sub1$Intermediate.Region.Name)
sub1$Intermediate.Region.Name<-revalue(sub1$Intermediate.Region.Name, c("Eastern Africa"="Eastern", "Middle Africa"="Middle", 
                                                                        "Southern Africa"="Southern", "Western Africa"="Western"))
sub1$Intermediate.Region.Name<-factor(sub1$Intermediate.Region.Name, levels=c("Northern","Western", "Middle", "Eastern", "Southern"))
sub1<-sub1[complete.cases(sub1$`Country of origin`),]

k3p<-ggplot(sub1, aes(x=sample, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  theme_cowplot() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(), 
        axis.title.y=element_blank(),
        strip.background = element_blank(), strip.text = element_blank(), 
        axis.title.x = element_blank(),
        axis.line=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  xlab("accession") +
  scale_fill_manual(values=pal) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(.~Region.Name+Intermediate.Region.Name, scales="free_x", space="free_x")

##########
## k = 4
k4<-df[df$K==4,]

## produce sort order
ord<-k4[,c("sample", "value", "variable")]
ord<-spread(ord, variable, value)
maxval<-apply(ord[,2:ncol(ord)], 1, max)
matchval<-vector(length=nrow(ord))
for(j in 1:nrow(ord)) matchval[j] <- match(maxval[j], ord[j, ])
ord$maxval<-maxval
ord$matchval<-matchval
ord<-ord[with(ord, order(matchval, -maxval)),]

k4$sample<-factor(k4$sample, levels=ord$sample)

## add regional data
k4<-merge(k4, samps, by="sample", all.x=T)

## split plots by region
## plot for Africa
sub1<-k4
sub1$Intermediate.Region.Name<-ifelse(sub1$Sub.region.Name=="Northern Africa", "Northern", sub1$Intermediate.Region.Name)
sub1$Intermediate.Region.Name<-revalue(sub1$Intermediate.Region.Name, c("Eastern Africa"="Eastern", "Middle Africa"="Middle", 
                                                                        "Southern Africa"="Southern", "Western Africa"="Western"))
sub1$Intermediate.Region.Name<-factor(sub1$Intermediate.Region.Name, levels=c("Northern","Western", "Middle", "Eastern", "Southern"))
sub1<-sub1[complete.cases(sub1$`Country of origin`),]

k4p<-ggplot(sub1, aes(x=sample, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  theme_cowplot() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(), 
        axis.title.y=element_blank(),
        strip.background = element_blank(), strip.text = element_blank(), 
        axis.title.x = element_blank(),
        axis.line=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  ylab("Proportion of ancestry") +
  xlab("accession") +
  scale_fill_manual(values=pal) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(.~Region.Name+Intermediate.Region.Name, scales="free_x", space="free_x")

## align pop structure panels
bp<-k2p/k3p/k4p 
print(bp)
###########
###########
## produce pca colored by cluster
# import data
## pc scores
df<-read.delim("../data/scores_pc.txt")
df$sample<-row.names(df)

## cowpea passport info
df1<-read_excel("../data/IITA-Cowpea collection.xls")
df1<-df1[,c("Accession name", "Latitude", "Longitude", "Biological status of accession")]
names(df1)[1]<-"sample"

## cluster assignments
df2<-read.delim("../data/cluster_assignments.txt")
df2<-as.data.frame(unique(df2[,c("sample", "K", "assignment")]))

## merge data frames for plotting
df<-merge(df, df1, by="sample", all.x=T)
df1<-merge(df, df2, by="sample", all.x=T)

df2<-read.delim("../data/prop_var_pc.txt")

## k = 2
df1$assignment<-factor(df1$assignment, levels=c("Cluster1", "Cluster2", "Cluster3", "Cluster4",
                                                "Admixed"))
k2<-df1[df1$K==2,]

k2p1<-ggplot() +
  geom_point(data=k2[k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2), color="gray69") +
  geom_point(data=k2[!k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2, color=assignment)) +
  geom_point() +
  theme_cowplot() +
  xlab("PC 1") +
  ylab("PC 2") +
  scale_color_manual(values=pal) +
  theme(legend.position="none", 
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("K = 2")
  
p1<-ggplot() +
  geom_point(data=k2[k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2), color="gray69") +
  geom_point(data=k2[!k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2, color=assignment)) +
  geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(df2[1,3]*100), "%)" )) +
  ylab(paste0("PC 2 (", round(df2[2,3]*100), "%)" )) +
  scale_color_manual(values=pal) +
  theme(legend.position="none")

ggsave("~/Desktop/pcak4.pdf", p1, height=4, width=4)


## k = 3
k2<-df1[df1$K==3,]

k3p1<-ggplot() +
  geom_point(data=k2[k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2), color="gray69") +
  geom_point(data=k2[!k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2, color=assignment)) +
  geom_point() +
  theme_cowplot() +
  xlab("PC 1") +
  ylab("PC 2") +
  scale_color_manual(values=pal) +
  theme(legend.position="none", 
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("K = 3")


## k = 4
k2<-df1[df1$K==4,]

k4p1<-ggplot() +
  geom_point(data=k2[k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2), color="gray69") +
  geom_point(data=k2[!k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2, color=assignment)) +
  geom_point() +
  theme_cowplot() +
  xlab("PC 1") +
  ylab("PC 2") +
  scale_color_manual(values=pal) +
  theme(legend.position="none", 
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("K = 4")


k2p1/k3p1/k4p1
##########
##########
# prop variance PCs 1-10
# import data
df<-read.delim("../data/prop_var_pc.txt")
df<-df[1:10,]

pv<-ggplot(df, aes(x=pc, y=proportion_variance*100)) + 
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("Principal component") +
  ylab("% variance explained") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0, 10,5), limits=c(0, 10))

#########
#########
# capture legend for plotting
legend <- cowplot::get_legend(k4p + theme(legend.position="right"))
leg<-as_ggplot(legend)
print(leg)

##########
##########
# put it all together
#row1<-(dK + k2p + k2p1) + plot_layout(widths=c(1,3,1))
#row2<-(pv + k3p + k3p1) + plot_layout(widths=c(1,3,1))
#row3<-(leg + k4p + k4p1) + plot_layout(widths=c(1,3,1))

#fig2<-row1/row2/row3 + plot_layout(heights=c(1,1,1))

#ggsave("~/Desktop/Fig2_DRAFT.pdf", fig2, height=8, width=12)
#####
# try layout
layout<-"
ABBBBC
DEEEEF
GHHHHI
"

fig2<-dK + k2p + k2p1 + pv + k3p + k3p1 + leg + k4p + k4p1 + plot_layout(design = layout)
ggsave("~/Desktop/Fig2_DRAFT.pdf", fig2, height=8, width=12)
