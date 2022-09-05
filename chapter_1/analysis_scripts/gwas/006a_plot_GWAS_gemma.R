#!/usr/bin/env Rscript
# Manhattan and QQ plots of GEMMA results

options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
args

library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table)

# set variables
FileName=args[1]

# set output file name 
OutName1<-unlist(strsplit(FileName, ".", fixed=T))[5]
OutName<-paste0("../..", OutName1)
  
# import data
df<-fread(FileName)

# Bonferroni threshold 
threshold <- 0.05/nrow(df)
  
# parse locus
names(df)[1]<-"CHR"
  
# following code adapted from:
# https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
# format for plotting
df$BP<-as.numeric(df$ps)
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
  
result<-result %>% filter(-log10(p_lrt)>2)

axisdf = result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Manhattan plot 
g<-ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
    
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=-log10(threshold)), color = "firebrick1", linetype="dashed", alpha=0.7) +
    scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5)) +     # remove space between plot area and x axis
    
    # Custom the theme:
    theme_classic() +
    theme(legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text=element_text(size=16)) + 
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle(OutName1)
 
OutName1<-paste0(OutName, "_manhattan.jpeg")
ggsave(OutName1, g, width=8, height=3, units="in")
  
 ## qq plot
print("qq")
obs<--log10(sort(df$p_score, decreasing=F))
exp<--log10(ppoints(length(obs)))
m<-max(obs) + 0.5
df2<-as.data.frame(cbind(obs, exp))
df2$sig<-ifelse(df2$obs > -log10(threshold), "yes", "no")
  
g<-ggplot(df2) + 
  geom_point(aes(x=exp, y=obs, color=sig), alpha=0.5) + 
  xlim(0, m) + 
  ylim(0, m) + 
  theme_classic() + theme(text=element_text(size=16)) + 
  ylab(expression(Observed~~-log[10](italic(p)))) +
  xlab(expression(Expected~~-log[10](italic(p)))) +
  theme(legend.position="none") + 
  scale_color_manual(values=c("black", "red")) +
  geom_abline(slope=1, intercept=0, color="blue", linetype="dashed", alpha=0.5)
  
OutName2<-paste0(OutName, "_qq.jpeg")
ggsave(OutName2, g, height = 4, width = 4)
