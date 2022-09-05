#!/usr/bin/env Rscript
# stacked manhattans
# cjfiscus
# 2022-06-15

setwd("~/Desktop")

library(pacman)
p_load(ggplot2, dplyr, data.table, tidyr, patchwork, cowplot, viridis, 
       ggrepel, stringr, Sushi, scales)

pheno_lst<-read.delim("MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/DATA/master_phenotype_table.txt")

lst<-c("14C", "15C", "16C", "18C", 
       "22C", "21C", "19C", "20C")

num<-1

marklst<-read.delim("tag_snps_focal.txt")
names(marklst)[1]<-"rs"
marklst$label<-marklst$rs

for (i in lst){
  # set output nm and import data
  nm<-pheno_lst[pheno_lst$ID==i,]$PHENOTYPE
  
  file<-paste0("MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/DATA/gwas/", i, ".assoc.txt")
  df<-fread(file)
  
  # set threshold
  threshold<-0.05/nrow(df)
  
  # prep for plotting
  names(df)[1]<-"CHR"
  df$BP<-as.numeric(df$ps)
  
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
  
  sub<-marklst[marklst$id==i,]
  result<-merge(result, sub, by="rs", all.x=T)
  
  # manhattan
  p<-ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
    
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5) +
    geom_point(data=result[!is.na(result$label),],color="firebrick1") +
    geom_hline(aes(yintercept=-log10(threshold)), color = "chocolate1", linetype="dashed", alpha=0.7) +
    scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0, 0)) +
    scale_y_continuous(breaks= pretty_breaks(n=5), expand=expansion(mult=c(0.05,0.10))) +
    #scale_y_continuous(expand = c(0, 1.1)) +     # remove space between plot area and x axis
    
    # Custom the theme:
    theme_cowplot() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle(nm)
  
  assign(paste0("p",num), p)
  num<-num+1
}
#########
# last manhattans with labels
i<-"97"
nm<-pheno_lst[pheno_lst$ID==i,]$PHENOTYPE

file<-paste0("MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/DATA/gwas/", i, ".assoc.txt")
df<-fread(file)

# set threshold
threshold<-0.05/nrow(df)

# prep for plotting
names(df)[1]<-"CHR"
df$BP<-as.numeric(df$ps)

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

sub<-marklst[marklst$id==i,]
result<-merge(result, sub, by="rs", all.x=T)

# manhattan
p10<-ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.5) +
  geom_point(data=result[!is.na(result$label),],color="firebrick1") +
  geom_hline(aes(yintercept=-log10(threshold)), color = "chocolate1", linetype="dashed", alpha=0.7) +
  scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
  
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
  ylab(expression(-log[10](italic(p)))) +
  ggtitle(nm) +
  xlab("Chromosome")
#########
i<-"24C"
nm<-pheno_lst[pheno_lst$ID==i,]$PHENOTYPE

file<-paste0("MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/DATA/gwas/", i, ".assoc.txt")
df<-fread(file)

# set threshold
threshold<-0.05/nrow(df)

# prep for plotting
names(df)[1]<-"CHR"
df$BP<-as.numeric(df$ps)

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

sub<-marklst[marklst$id==i,]
result<-merge(result, sub, by="rs", all.x=T)

# manhattan
p11<-ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.5) +
  geom_point(data=result[!is.na(result$label),],color="firebrick1") +
  geom_hline(aes(yintercept=-log10(threshold)), color = "chocolate1", linetype="dashed", alpha=0.7) +
  scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
  
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
  ylab(expression(-log[10](italic(p)))) +
  ggtitle(nm) +
  xlab("Chromosome")
#########


# stack manhattans
col1<-p1 / p2 / p3 / p4 / p10
col2<-p5/p6/p7/p8/p11

ggsave("~/Desktop/col1.png", col1, height=6, width=6, dpi=300)
ggsave("~/Desktop/col1.pdf", col1, height=6, width=6)
ggsave("~/Desktop/col2.png", col2, height=6, width=6, dpi=300)
ggsave("~/Desktop/col2.pdf", col2, height=6, width=6)
##########
##########

## plot for zoomed manhattan plot
df1<-read.delim("profiled_loci.txt")

for (i in 1:nrow(df1)){
  sub<-df1[i,]
  chrom<-sub$chrom
  chromstart<-sub$chromstart-1000
  chromend<-sub$chromend+1000
  #chromend<-ifelse(chromend-bp > 500000, bp+500000, chromend)
  #chromstart<-ifelse(bp-chromstart > 500000, bp-500000, chromstart)
  tag<-sub$tag
  bp<-sub$bp
  cand<-unlist(str_split(sub$cand, ","))
  file<-paste0("MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/DATA/gwas/", sub$i, ".assoc.txt")
  
  df<-fread(file)
  threshold<-0.05/nrow(df)
  
  ### filter gwas snps
  df<-df[df$chr==chrom,]
  df<-df[df$ps >= chromstart,]
  df<-df[df$ps<=chromend,]
  
  ## import ld and filter
  ld<-fread("tag_snps.ld")
  ld<-ld[ld$SNP_A == tag,]
  ld<-ld[,c("SNP_B", "R2")]
  names(ld)[1]<-"rs"
  df<-merge(df, ld, by="rs", all.x=T)
  df$ltag<-ifelse(df$rs==tag, "y", "n")
  
  # plot for ld in region
  pa1<-ggplot() + 
    geom_point(data=df[df$ltag=="n",], aes(x=ps, y=-log10(p_lrt), color=R2)) +
    geom_point(data=df[df$ltag=="y",], aes(x=ps, y=-log10(p_lrt)), size=2, 
               shape=17, color="#7A0403FF") +
    geom_hline(aes(yintercept=-log10(threshold)), color="chocolate1", linetype="dashed") +
    geom_vline(aes(xintercept=bp), color="dimgrey", linetype="dashed") +
    theme_cowplot() +
    scale_color_viridis(option="H", limits = c(0.2, 1)) +
    ylab(expression(-log[10](italic(p)))) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + ggtitle(tag) +
    xlim(c(chromstart, chromend)) #+
    scale_x_continuous(labels=function(x)x/1000000)
  
  # plot for gene track
  genes<-read.table("genes.bed")
  names(genes)<-c("chrom", "start", "stop", "gene", "score", "strand")
  
  ## subset genes in region
  sub<-genes[genes$chrom==chrom,]
  sub<-sub[sub$start >= chromstart,]
  sub<-sub[sub$stop <= chromend,]
  
  ## define coords for arrows
  sub$astart<-ifelse(sub$strand=="+", sub$start, sub$stop)
  sub$astop<-ifelse(sub$strand=="+", sub$stop, sub$start)
  sub$labx<-((sub$stop-sub$start)/2) + sub$start
  sub$rank<-1:nrow(sub)
  sub$iscand<-ifelse(sub$gene %in% cand, "y", "n")
  
  pa2<-ggplot(sub) + 
    geom_vline(aes(xintercept=bp), color="dimgrey", linetype="dashed") +
    geom_segment(data=sub[sub$iscand=="n",],aes(x=astart, xend=astop, y=rank, yend=rank), 
                 arrow = arrow(length = unit(0.15, "cm"), type="closed"), size=1.2) + 
    geom_segment(data=sub[sub$iscand=="y",],aes(x=astart, xend=astop, y=rank, yend=rank), 
                 arrow = arrow(length = unit(0.15, "cm"), type="closed"), size=1.2, color="firebrick1") +
    geom_text(data=sub[sub$iscand=="y",], aes(x=labx, y=(rank + max(sub$rank)/4), label=gene, fontface = "italic")) +
    theme_cowplot() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.line.y=element_blank()) +
    xlab(paste0(chrom," position (Mbp)")) + xlim(c(chromstart, chromend)) +
    scale_y_continuous(expand=expansion(mult=c(0.20,0.20))) #+
    scale_x_continuous(labels=function(x)x/1000000)
  
  plot<-(pa1 + theme(legend.position="none"))/pa2 + plot_layout(heights=c(3,1.5))
  ggsave(paste0(i, "_locus.pdf"), plot, height=3, width=3)
}

## export legend
leg<-get_legend(pa1)
ggsave("legend.pdf", leg, height=2, width=2)
