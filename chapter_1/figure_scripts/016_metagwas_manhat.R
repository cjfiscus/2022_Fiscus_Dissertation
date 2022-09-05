#!/usr/bin/env Rscript
# meta-GWAS manhattans
# cjfiscus
# 2021-12-09

# load libs
library(pacman)
p_load(data.table, ggplot2, tidyr, dplyr, ggpubr, patchwork, cowplot, ggupset, 
       tidyverse, dplyr, ggnewscale)

# meta-GWA scans (A-D)
## A Genes
### shared meta-GWA hits to highlight
df1<-read.delim("~/Desktop/MANUSCRIPTS/athal_data/metagwas_top1326_unique_count.txt")
df1<-df1[df1$Freq > 1,]

##########
### scan for retrotransposons
file="~/Desktop/MANUSCRIPTS/athal_data/retrotrans_meta.txt"
k=88

# parse variables
df<-fread(file)
prefix<-unlist(strsplit(file, "_"))[1]

# calculate lambda
chisq<-median(df$chisq)/qchisq(0.5, k)

# correct chisq values
df$chisq_adj<-df$chisq/chisq

# calculate threshold
df<-df[order(df$chisq_adj, decreasing=T),]
df$rank<-1:nrow(df)
df$highlight<-ifelse(df$rank <= 1326, "yes", "no")
df<-merge(df, df1, by="snp_id", all.x=T)
df$Freq<-ifelse(df$highlight=="yes" & is.na(df$Freq), 1, df$Freq)

# qqplots before and after genomic inflation correction
param<-list(df=2*k)

# plot manhattan 
df<-df %>% separate(snp_id, c("CHR", "BP"))
df$CHR<-factor(df$CHR)
df$BP<-as.numeric(as.character(df$BP))

result <- df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

result<-result %>% filter(chisq_adj > quantile(df$chisq_adj, probs=c(0.90)))
axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
result$SNPID<-paste0(result$CHR, ":", result$BP)

### manhattan 
mypal<-c("firebrick1", "dodgerblue1")
names(mypal)<-c("3", "2")

g1<-ggplot() +
  
  # Show all points
  geom_point(data=subset(result, highlight=="no"), aes(x=BPcum, y=chisq_adj, color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c("gray19", "gray48"), 22 )) +
  
  # snp highlight
  new_scale_color() +
  geom_point(data=subset(result, highlight=="yes" & Freq > 1), aes(x=BPcum, y=chisq_adj, color=as.factor(Freq)), size=3, alpha=0.5) +
  
  # significance line
  #geom_hline(aes(yintercept=threshold3), color="dimgrey", linetype="dashed") +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand=expansion(mult=c(0.01,0.01))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.01))) +
  
  # Custom the theme:
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("") +
  ylab(expression(italic(chi^2))) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values=mypal) +
  ggtitle("Retrotransposon")
#ggsave("retrotrans.tiff", g1, height=1.5, width=10)

###########
### scan for DNA transposons
### scan for retrotransposons
file="~/Desktop/MANUSCRIPTS/athal_data/dnatrans_meta.txt"
k=31

# parse variables
df<-fread(file)
prefix<-unlist(strsplit(file, "_"))[1]

# calculate lambda
chisq<-median(df$chisq)/qchisq(0.5, k)

# correct chisq values
df$chisq_adj<-df$chisq/chisq

# calculate threshold
df<-df[order(df$chisq_adj, decreasing=T),]
df$rank<-1:nrow(df)
df$highlight<-ifelse(df$rank <= 1326, "yes", "no")
df<-merge(df, df1, by="snp_id", all.x=T)
df$Freq<-ifelse(df$highlight=="yes" & is.na(df$Freq), 1, df$Freq)

# qqplots before and after genomic inflation correction
param<-list(df=2*k)

# plot manhattan 
df<-df %>% separate(snp_id, c("CHR", "BP"))
df$CHR<-factor(df$CHR)
df$BP<-as.numeric(as.character(df$BP))

result <- df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

result<-result %>% filter(chisq_adj > quantile(df$chisq_adj, probs=c(0.90)))
axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
result$SNPID<-paste0(result$CHR, ":", result$BP)

### manhattan 
mypal<-c("firebrick1", "darkorange")
names(mypal)<-c("3", "2")

mypal2<-c(5, 2)
names(mypal2)<-c("3", "2")

g2<-ggplot() +
  
  # Show all points
  geom_point(data=subset(result, highlight=="no"), aes(x=BPcum, y=chisq_adj, color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c("gray19", "gray48"), 22 )) +
  
  # snp highlight
  new_scale_color() +
  geom_point(data=subset(result, highlight=="yes" & Freq > 1), aes(x=BPcum, y=chisq_adj, color=as.factor(Freq)), size=3, alpha=0.5) +
  
  # significance line
  #geom_hline(aes(yintercept=threshold3), color="dimgrey", linetype="dashed") +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand=expansion(mult=c(0.01,0.01))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.01))) +
  
  # Custom the theme:
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("") +
  ylab(expression(italic(chi^2))) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values=mypal)
#ggsave("dnatrans.jpeg", g2, height=1.5, width=10)

#########
### scan for simple repeats
### scan for retrotransposons
file="~/Desktop/MANUSCRIPTS/athal_data/sr_meta.txt"
k=60

# parse variables
df<-fread(file)
prefix<-unlist(strsplit(file, "_"))[1]

# calculate lambda
chisq<-median(df$chisq)/qchisq(0.5, k)

# correct chisq values
df$chisq_adj<-df$chisq/chisq

# calculate threshold
df<-df[order(df$chisq_adj, decreasing=T),]
df$rank<-1:nrow(df)
df$highlight<-ifelse(df$rank <= 1326, "yes", "no")
df<-merge(df, df1, by="snp_id", all.x=T)
df$Freq<-ifelse(df$highlight=="yes" & is.na(df$Freq), 1, df$Freq)

# qqplots before and after genomic inflation correction
param<-list(df=2*k)

# plot manhattan 
df<-df %>% separate(snp_id, c("CHR", "BP"))
df$CHR<-factor(df$CHR)
df$BP<-as.numeric(as.character(df$BP))

result <- df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

result<-result %>% filter(chisq_adj > quantile(df$chisq_adj, probs=c(0.90)))
axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
result$SNPID<-paste0(result$CHR, ":", result$BP)

### manhattan 
mypal<-c("firebrick1", "darkorange")
names(mypal)<-c("3", "2")

mypal2<-c(5, 2)
names(mypal2)<-c("3", "2")

g3<-ggplot() +
  
  # Show all points
  geom_point(data=subset(result, highlight=="no"), aes(x=BPcum, y=chisq_adj, color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c("gray19", "gray48"), 22 )) +
  
  # snp highlight
  new_scale_color() +
  geom_point(data=subset(result, highlight=="yes" & Freq > 1), aes(x=BPcum, y=chisq_adj, color=as.factor(Freq)), size=3, alpha=0.5) +
  
  # significance line
  #geom_hline(aes(yintercept=threshold3), color="dimgrey", linetype="dashed") +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand=expansion(mult=c(0.01,0.01))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.01))) +
  
  # Custom the theme:
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("") +
  ylab(expression(italic(chi^2))) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values=mypal)
#ggsave("sr.tiff", g3, height=1.5, width=10)
########

# panels A-D
plot<-g1/g2/g3
ggsave("~/Desktop/meta.png", plot, height=6, width=12)
#ggsave("16_metagwas_manhat.pdf", plot, height=6, width=12)