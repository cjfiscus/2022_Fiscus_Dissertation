#!/usr/bin/env Rscript
# Figure 4
# cjfiscus
# 2022-07-26

library(pacman)
p_load(ggplot2, cowplot, ggpubr, data.table, dplyr, scales, ggrepel, tidyverse,
       patchwork)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/SCRIPTS/")

# A bioclim PCA 
df<-read.delim("../data/bioclim_pc_scores.txt")
df1<-read.delim("../data/cluster_assignments.txt")
df1<-df1[df1$K==2,]
df1<-as.data.frame(unique(df1[,c("sample", "assignment")]))
df<-merge(df, df1, by="sample")                   
pervar<-read.delim("../data/bioclim_per_var_exp.txt")

pal<-c("#4DBBD5FF", "#E64B35FF", "dimgrey")
names(pal)<-c("Cluster1", "Cluster2", "Admixed")

p1<-ggplot(df, aes(x=PC1, y=PC2, color=assignment)) + 
  geom_point(alpha=0.5) + theme_cowplot() +
  scale_color_manual(values=pal) +
  xlab(paste0("PC 1 (", round(pervar[1,2]), "%)")) +
  ylab(paste0("PC 2 (", round(pervar[2,2]), "%)"))

# B PC loadings
df<-read.delim("../data/bioclim_pc_loadings.txt")
pal2<-c("aquamarine4", "chocolate1")
names(pal2)<-c("Precipitation", "Temperature")
p2<-ggplot(df, aes(x=PC1, y=PC2, label=lab, color=type)) + 
  geom_text_repel() + 
  geom_point() + theme_cowplot() +
  scale_color_manual(values=pal2)

###########
# PC 1 vs flowering time
df<-read.delim("../data/rad_data_samples_cw.txt")
names(df)[1]<-"sample"
df1<-read.delim("../data/bioclim_pc_scores.txt")
df<-merge(df, df1, by="sample")
df1<-read.delim("../data/cluster_assignments.txt")
df1<-df1[df1$K==2,]
df1<-as.data.frame(unique(df1[,c("sample", "assignment")]))
df<-merge(df, df1, by="sample")

lst<-c(13, 16, 19, 22, 23, 24)

cs<-data.frame()
for (i in lst){
  nm<-names(df)[i]
  c<-cor.test(df[,60], as.numeric(df[,i]), use="pairwise.complete.obs")
  c2<-cor.test(df[,61], as.numeric(df[,i]), use="pairwise.complete.obs")
  
  cs<-as.data.frame(rbind(cs, cbind("pc1", nm, c$p.value, c$estimate),
  cbind("pc2", nm, c2$p.value, c2$estimate)))
}

p3<-ggplot(df, 
       aes(x=PC1, y=FlowStart.day_after_sowing._Fam_mean)) + 
  geom_point(aes(color=assignment), alpha=0.5) +
  geom_smooth(method="lm") + theme_cowplot() +
  scale_color_manual(values=pal) + ylab("Days to flowering") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x.npc = "centre")
##########
# PC 1 GWAS
df<-fread("../data/bioclim_PC1.assoc.txt")

# set threshold
threshold<-0.05/nrow(df)

# prep for plotting
names(df)[1]<-"CHR"
df$BP<-as.numeric(df$ps)
df$CHR<-gsub("SCF_", "", df$CHR)
df$CHR<-factor(df$CHR, levels=1:16)

result <- df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

lst<-read.delim("../data/highlight_snps.txt")
lst<-lst$x

#lst<-c("SCF_14:15377114", "SCF_9:8559318", "SCF_14:15183887")

result$rs<-paste0("SCF_", result$CHR, ":", result$ps)
result$lab<-ifelse(result$rs %in% lst, "y", "n")

# manhattan
p4<-ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.5) +
  geom_point(data=result[result$lab=="y",],color="firebrick1") +
  geom_hline(aes(yintercept=-log10(threshold)), 
             color = "chocolate1", linetype="dashed", alpha=0.7) +
  scale_color_manual(values = c(rep(c("#7570B3", "#ABA4EB"), 4), 
                                rep(c("#1B9E77", "#60D6AB"), 4))) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0, 0)) +
  scale_y_continuous(breaks= pretty_breaks(n=5), expand=expansion(mult=c(0.05,0.10))) +
  #scale_y_continuous(expand = c(0, 1.1)) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab(expression(-log[10](italic(p)))) + xlab("chromosome")

###########
# PC 2 GWAS
df<-fread("../data/bioclim_PC2.assoc.txt")

# set threshold
threshold<-0.05/nrow(df)

# prep for plotting
names(df)[1]<-"CHR"
df$BP<-as.numeric(df$ps)
df$CHR<-gsub("SCF_", "", df$CHR)
df$CHR<-factor(df$CHR, levels=1:16)

result <- df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

result$rs<-paste0("SCF_", result$CHR, ":", result$ps)
result$lab<-ifelse(result$rs %in% lst, "y", "n")

# manhattan
p5<-ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.5) +
  geom_point(data=result[result$lab=="y",],color="firebrick1") +
  geom_hline(aes(yintercept=-log10(threshold)), 
             color = "chocolate1", linetype="dashed", alpha=0.7) +
  scale_color_manual(values = c(rep(c("#7570B3", "#ABA4EB"), 4), 
                                rep(c("#1B9E77", "#60D6AB"), 4))) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0, 0)) +
  scale_y_continuous(breaks= pretty_breaks(n=5), expand=expansion(mult=c(0.05,0.10))) +
  #scale_y_continuous(expand = c(0, 1.1)) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab(expression(-log[10](italic(p)))) + xlab("chromosome")

##########
# PC1 go
df<-read.delim("../data/GWAS_bioclimpc1_closest_GO.txt")

df$label<-gsub('(.{1,60})(\\s|$)', '\\1\n', df$Term)

p6a<-ggplot(df, aes(x=reorder(label, enrich, FUN="identity"), 
               y=enrich)) + geom_bar(stat="identity") + 
  theme_cowplot() + ylab("Enrichment") + coord_flip() + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())
p6b<-ggplot(df, aes(x=reorder(label, enrich, FUN="identity"), 
                    y=enrich)) + 
  theme_cowplot() + ylab("Enrichment") + coord_flip() + 
  theme(axis.title.y=element_blank())

##########
# PC1 go
df<-read.delim("../data/GWAS_bioclimpc2_closest_GO.txt")
df<-df[1:20,]

df$label<-gsub('(.{1,60})(\\s|$)', '\\1\n', df$Term)

p7a<-ggplot(df, aes(x=reorder(label, enrich, FUN="identity"), 
                    y=enrich)) + geom_bar(stat="identity") + 
  theme_cowplot() + ylab("Enrichment") + coord_flip() + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())
p7b<-ggplot(df, aes(x=reorder(label, enrich, FUN="identity"), 
                    y=enrich)) + 
  theme_cowplot() + ylab("Enrichment") + coord_flip() + 
  theme(axis.title.y=element_blank())
##########
# export panels
ggsave("~/Desktop/A.pdf", p1 + theme(legend.position="none"), height=3, width=3)
ggsave("~/Desktop/B.pdf", p2 + theme(legend.position="none"), height=3, width=3)
ggsave("~/Desktop/C.pdf", p3 + theme(legend.position="none"), height=3, width=3)
ggsave("~/Desktop/D.pdf", p4 + theme(legend.position="none"), height=2, width=10)
ggsave("~/Desktop/E.pdf", p5 + theme(legend.position="none"), height=2, width=10)
ggsave("~/Desktop/F1.pdf", p6a + theme(legend.position="none"), height=3, width=1.5)
ggsave("~/Desktop/F2.pdf", p6b + theme(legend.position="none"), height=3, width=6)
ggsave("~/Desktop/G1.pdf", p7a + theme(legend.position="none"), height=3, width=1.5)
ggsave("~/Desktop/G2.pdf", p7b + theme(legend.position="none"), height=3, width=6)
leg<-get_legend(p1)
ggsave("~/Desktop/leg1.pdf", leg, height=2, width=2)
leg<-get_legend(p2)
ggsave("~/Desktop/leg2.pdf", leg, height=2, width=2)
########## OLD CODE
##########
# gene plot for tag snp
df<-fread("../data/bioclim_PC1.assoc.txt")
df<-as.data.frame(df)

chrom<-"SCF_14"
bp<-15183887
chromstart<-bp-10000
chromend<-bp+10000
cand<-"g39808"

# set threshold
threshold<-0.05/nrow(df)

## filter gwas snps
df<-df[df$chr==chrom,]
df<-df[df$ps >= chromstart,]
df<-df[df$ps <=chromend,]

## filter genes
genes<-read.table("../data/Cbp_genes.bed")
names(genes)<-c("chr", "start", "end", "strand", "id")
sub<-genes[genes$chr==chrom,]
sub<-sub[sub$start >=chromstart,]
sub<-sub[sub$end<=chromend,]
sub$astart<-ifelse(sub$strand=="+", sub$start, sub$end)
sub$astop<-ifelse(sub$strand=="+", sub$end, sub$start)
sub$labx<-((sub$end-sub$start)/2) + sub$start
sub$iscand<-ifelse(sub$id %in% cand, "y", "n")
sub$rank<-1:nrow(sub)

p5<-ggplot(sub) + 
  geom_vline(aes(xintercept=bp), color="dimgrey", linetype="dashed") +
  geom_segment(data=sub[sub$iscand=="n",], 
                aes(x=astart, xend=astop, y=rank, yend=rank),
               arrow = arrow(length = unit(0.15, "cm"), type="closed"), 
               size=1.2) +
  geom_segment(data=sub[sub$iscand=="y",], 
               aes(x=astart, xend=astop, y=rank, yend=rank),
               arrow = arrow(length = unit(0.15, "cm"), type="closed"), 
               size=1.2,color="firebrick1") +
  geom_text(data=sub, aes(x=labx, y=(rank + max(rank)/4), label=id, fontface = "italic")) +
  theme_cowplot() + xlim(chromstart, chromend)

p6<-ggplot(df, aes(x=ps, y=-log10(p_lrt))) + geom_point() +
  theme_cowplot() + xlim(chromstart,chromend) + 
  geom_hline(aes(yintercept=-log10(threshold)), color = "chocolate1", linetype="dashed", alpha=0.7) +
  ylab(expression(-log[10](italic(p)))) +
  geom_vline(aes(xintercept=bp), color="dimgrey", linetype="dashed")
##########

df<-read.delim("~/Desktop/2022-07-25_CBP_GWAS_REGIONS - Sheet1.tsv")
df1<-df[df$ID=="bioclim_PC2",]
df2<-df1[nchar(df1$NOTE) > 0,]
write.table(df2, "~/Desktop/pc2_cands.txt", sep="\t", quote=F, row.names=F)

df<-read.delim("~/Desktop/post_gwas/GWAS_bioclimpc1_closest_GO.txt")
df1<-read.delim("~/Desktop/post_gwas/GWAS_bioclimpc2_closest_GO.txt")

table(df1$GO.ID %in% df$GO.ID)
df2<-df1[!df1$GO.ID %in% df$GO.ID,]
