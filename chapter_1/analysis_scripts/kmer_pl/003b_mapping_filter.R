#!/usr/bin/env Rscript
# Contamination filter
# cjfiscus

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2)

args = commandArgs(trailingOnly=TRUE)

# Identify libraries with < 95% reads mapping
# Produce plots showing effects of mapping

# read in data
df<-read.table(args[1])
df$V2<-as.numeric(as.character(df$V2))
df<-na.omit(df)

# distribution plot
g<-ggplot(df) + geom_density(aes(x=V2)) + xlab("% Reads Mapped") + 
	theme_classic() + theme(text=element_text(size=18))
ggsave("../results/filters/ref_genome_map.png", g, width=5, height=5, units="in")

# filter at 90% 
low<-90
df1<-df[df$V2<low,]
write.table(df1$V1, args[2], sep="\t", quote=F, row.names=F, col.names = F)

g<-ggplot(df) + geom_density(aes(x=V2)) + xlab("% Reads Mapped") + 
	theme_classic() + theme(text=element_text(size=18)) + 
	geom_vline(aes(xintercept=low), color="red", linetype="dashed", alpha=0.5)
ggsave("../results/filters/ref_genome_map_cutoff.png", g, width=5, height=5, units="in")

# data post-filtering 
g<-ggplot(df) + geom_density(aes(x=V2)) + xlab("% Reads Mapped") + 
	theme_classic() + theme(text=element_text(size=18)) + 
	scale_x_continuous(limits = c(low, 100))
ggsave("../results/filters/ref_genome_map_limit.png", g, width=5, height=5, units="in")
