#!/usr/bin/env Rscript

library(pacman)
p_load(cowplot, GGally)

# read in and format data
df<-read.delim("../results/table_12mers.txt")
names(df)[names(df) == "X1741_Norm"]<-"Norm"
names(df)[names(df) == "X1741_Raw"]<-"Raw"
names(df)[names(df) == "Illumina.MiSeq"]<-"MiSeq"

# write out plot
g<-ggcorr(df[,2:ncol(df)], method = c("all.obs", "spearman"), 
       label=TRUE, label_round = 2, size = 3, layout.exp = 1, hjust=0.75, angle=0)
save_plot("../results/corr.jpeg", g)
