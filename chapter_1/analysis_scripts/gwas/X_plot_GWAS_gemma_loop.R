#!/usr/bin/env Rscript
# Manhattan and QQ plots of GEMMA results

options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
args

library(pacman)
p_load(ggplot2, dplyr, tidyr, data.table, gridExtra)

addPlot<-function(FileName){
  # set output file name 
  OutName1<-unlist(strsplit(FileName, ".", fixed=T))[5]
  OutName1<-unlist(strsplit(OutName1, "/"))[5]

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
  ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
    
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
  
  #print(g)
}

# set variables
df<-read.table(args[1])
lst<-df$V1

test<-lapply(lst, addPlot)
ggsave(args[2], marrangeGrob(grobs=test, nrow=4, ncol=1))
