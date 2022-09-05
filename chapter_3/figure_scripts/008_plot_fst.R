#!/usr/bin/env Rscript
# plot fst
# cjfiscus
# 2022-03-07

setwd("~/Desktop")

library(pacman)
p_load(ggplot2, cowplot, reshape2, tidyr, dplyr)

# import data
df<-read.delim("fst_plot/cluster1_vs_cluster2.windowed.weir.fst")
df$Cluster1_vs_Cluster2<-df$MEAN_FST
df<-df[,c("CHROM", "BIN_START", "BIN_END", "Cluster1_vs_Cluster2")]

df1<-read.delim("fst_plot/cluster1_vs_cluster3.windowed.weir.fst")
df1$Cluster1_vs_Cluster3<-df1$MEAN_FST
df1<-df1[,c("CHROM", "BIN_START", "BIN_END", "Cluster1_vs_Cluster3")]

df2<-read.delim("fst_plot/cluster2_vs_cluster3.windowed.weir.fst")
df2$Cluster2_vs_Cluster3<-df2$MEAN_FST
df2<-df2[,c("CHROM", "BIN_START", "BIN_END", "Cluster2_vs_Cluster3")]

# merge
m<-merge(df, df1, by=c("CHROM", "BIN_START", "BIN_END"))
m<-merge(m, df2, by=c("CHROM", "BIN_START", "BIN_END"))

m<-melt(m, id.vars=c("CHROM", "BIN_START", "BIN_END"))

## 1 vs 2. 

m1<-m[m$variable=="Cluster1_vs_Cluster2",]

names(m1)[1]<-"CHR"
m1$BP<-((m1$BIN_END-m1$BIN_START-1)/2)+m1$BIN_START

result <- m1 %>%
  group_by(CHR) %>%
  summarize(chr_len=max(BP)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(m1, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)
 axisdf=result %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2)

result$value<-ifelse(result$value < 0, 0, result$value)

winm<-mean(result$value)

p1<-ggplot(result, aes(x=BPcum, y=value)) +
   
   # Show all points
   geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
   scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
  geom_hline(aes(yintercept=winm), color = "firebrick1", linetype="solid", alpha=0.7) +
   # custom X axis:
   scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
   scale_y_continuous(breaks=seq(0,1, 0.25), limits=c(0,1)) +     # remove space between plot area and x axis
   
   # Custom the theme:
   theme_classic() +
   theme(legend.position="none",
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         text=element_text(size=16)) + 
   xlab("Chromosome") +
   ylab("Fst") +
  ggtitle("Cluster 1 vs. Cluster 2")
ggsave("008A_fst_cluster1_vs_cluster2.pdf", p1, height=4, width=12) 

##########
m1<-m[m$variable=="Cluster1_vs_Cluster3",]

names(m1)[1]<-"CHR"
m1$BP<-((m1$BIN_END-m1$BIN_START-1)/2)+m1$BIN_START

result <- m1 %>%
  group_by(CHR) %>%
  summarize(chr_len=max(BP)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(m1, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)
axisdf=result %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2)

result$value<-ifelse(result$value < 0, 0, result$value)
winm<-mean(result$value)
p1<-ggplot(result, aes(x=BPcum, y=value)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
  geom_hline(aes(yintercept=winm), color = "firebrick1", linetype="solid", alpha=0.7) +
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(breaks=seq(0,1, 0.25), limits=c(0,1)) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_classic() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(size=16)) + 
  xlab("Chromosome") +
  ylab("Fst") +
  ggtitle("Cluster 1 vs. Cluster 3")
ggsave("008B_fst_cluster1_vs_cluster3.pdf", p1, height=4, width=12) 
##########
# cluster 2 vs. 3
m1<-m[m$variable=="Cluster2_vs_Cluster3",]

names(m1)[1]<-"CHR"
m1$BP<-((m1$BIN_END-m1$BIN_START-1)/2)+m1$BIN_START

result <- m1 %>%
  group_by(CHR) %>%
  summarize(chr_len=max(BP)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(m1, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)
axisdf=result %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2)

result$value<-ifelse(result$value < 0, 0, result$value)
winm<-mean(result$value)
p1<-ggplot(result, aes(x=BPcum, y=value)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
  geom_hline(aes(yintercept=winm), color = "firebrick1", linetype="solid", alpha=0.7) +
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(breaks=seq(0,1, 0.25), limits=c(0,1)) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_classic() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(size=16)) + 
  xlab("Chromosome") +
  ylab("Fst") +
  ggtitle("Cluster 2 vs. Cluster 3")
ggsave("008C_fst_cluster2_vs_cluster3.pdf", p1, height=4, width=12) 
