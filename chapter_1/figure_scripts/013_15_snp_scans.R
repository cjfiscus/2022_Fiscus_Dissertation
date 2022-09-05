#!/usr/bin/env Rscript
# windowed SNPs
# cjfiscus
# 2021-12-09

library(pacman)
p_load(data.table, IRanges, GenomicRanges, ggplot2, cowplot, patchwork)

## load data
df<-fread("~/Desktop/MANUSCRIPTS/athal_data/snps_lowp.txt.gz")
tests<-1325631
alpha<-0.05/tests

## subset to sig snps
df<-df[,c("V1", "V2", "V4", "V15")]
names(df)<-c("ID", "CHR", "BP", "P")
df<-df[df$P < alpha,]

## remove old data
df<-df[!df$ID=="dist_lyrata",]
df<-df[!df$ID=="dist_rubella",]

rpt<-read.delim("~/Desktop/MANUSCRIPTS/athal_data/A_thaliana_Repeat_Database.txt")
df<-df[df$ID %in% rpt$ID,]

## save SNPs to GRanges
snps<-GRanges(seqnames = df$CHR, 
              ranges=IRanges(start=as.numeric(df$BP), end=as.numeric(df$BP)+1), id=df$ID)

## set seq lens
ranges<-read.delim("~/Desktop/MANUSCRIPTS/athal_data/TAIR10_ranges.txt")
ranges<-ranges[1:5,]
rng<-ranges$END
names(rng)<-ranges$CHR
seqlengths(snps)<-rng

## count hits per window
### define window
win<-100000
bins<-tileGenome(seqlengths(snps), tilewidth=win, cut.last.tile.in.chrom=TRUE)

### tally hits per window
hits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(snps))))
names(hits)[6]<-"snps"
hits$pos<-hits$start+(hits$end-hits$start)/2

feats<-read.delim("~/Desktop/MANUSCRIPTS/athal_data/features.txt")
names(feats)[1]<-"seqnames"

### plot
p1<-ggplot() +
  geom_rect(data=feats[feats$feat=="CEN",], 
            aes(xmin=start, xmax=end, ymin=0, ymax=Inf), fill="lightpink", 
            alpha=0.5) +
  geom_bar(data=hits, aes(x=pos, y=snps, group=seqnames), 
                        stat="identity", color="gray37") +
  facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
  xlab("window") + ylab("sig. SNPs") +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank()) +
  xlab("position (Mbp)")
ggsave("~/Desktop/13_gwas_hits_100000_ALL.pdf", p1, height=3, width=8)
##########

# plot hits per class
for (i in unique(rpt$CLASS)){
  sub<-rpt[rpt$CLASS==i,]
  snps_sub<-snps[snps$id %in% sub$ID,]
  
  subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(snps_sub))))
  names(subhits)[6]<-"snps"
  subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 
  
  p1<-ggplot() + 
    geom_rect(data=feats[feats$feat=="CEN",], 
              aes(xmin=start, xmax=end, ymin=0, ymax=Inf), fill="lightpink", 
              alpha=0.5) +
    geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
                      stat="identity", color="gray37") +
    facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
    xlab("window") + ylab("sig. SNPs") +
    scale_x_continuous(labels=function(x)x/1000000) +
    theme(
      strip.background = element_blank()) +
    xlab("position (Mbp)")
  assign(paste0("b", i), p1)
  ggsave(paste0("~/Desktop/13_gwas_hits_", 
                as.character(as.integer(win)), "_", i, ".pdf"), 
         p1, height=3, width=8)
}

out<-bRetrotransposon / `bDNA transposon` / bSatellite / `bSimple Repeat`
ggsave("~/Desktop/out2.pdf", out, height=10, width=10)

## produce windowed plots for subclasses
lst<-unique(rpt$SUBCLASS)
lst<-lst[lst != "Other"]
for (i in unique(lst)){
  sub<-rpt[rpt$SUBCLASS==i,]
  sub<-sub[complete.cases(sub[,1:4]),]
  snps_sub<-snps[snps$id %in% sub$ID,]
  
  subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(snps_sub))))
  names(subhits)[6]<-"snps"
  subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 
  
  p1<-ggplot() + geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
                          stat="identity", color="gray37") +
    facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
    xlab("window") + ylab("sig. SNPs") +
    scale_x_continuous(labels=function(x)x/1000000) +
    theme(
      strip.background = element_blank()) +
    xlab("position (Mbp)")
  ggsave(paste0("~/Desktop/15_gwas_hits_", 
                as.character(as.integer(win)), "_",unique(sub$CLASS), "_", i, ".pdf"), 
         p1, height=3, width=8)
}

# plots by transposon family 
lst<-unique(rpt$FAMILY)
lst<-lst[lst != "Other"]
for (i in unique(lst)){
  sub<-rpt[rpt$FAMILY==i,]
  sub<-sub[complete.cases(sub[,1:4]),]
  snps_sub<-snps[snps$id %in% sub$ID,]
  
  subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(snps_sub))))
  names(subhits)[6]<-"snps"
  subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 
  
  p1<-ggplot() + geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
                          stat="identity", color="gray37") +
    facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
    xlab("window") + ylab("sig. SNPs") +
    scale_x_continuous(labels=function(x)x/1000000) +
    theme(
      strip.background = element_blank()) +
    xlab("position (Mbp)")
  ggsave(paste0("~/Desktop/15_gwas_hits_", 
                as.character(as.integer(win)), "_", unique(sub$CLASS), "_", sub("/", "_", i), ".pdf"), 
         p1, height=3, width=8)
}

## Produce other plots
### satellites other
sub<-rpt[rpt$CLASS=="Satellite" & rpt$SUBCLASS=="Other",]
snps_sub<-snps[snps$id %in% sub$ID,]

subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(snps_sub))))
names(subhits)[6]<-"snps"
subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 

p1<-ggplot() + geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
                    stat="identity", color="gray37") +
  facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
  xlab("window") + ylab("sig. SNPs") +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank()) +
  xlab("position (Mbp)")
ggsave(paste0("~/Desktop/15_gwas_hits_", as.character(as.integer(win)), "Satellite_Other.pdf"),
       p1, height=3, width=8)

### DNA transposon Other
sub<-rpt[rpt$CLASS=="DNA transposon" & rpt$FAMILY=="Other",]
snps_sub<-snps[snps$id %in% sub$ID,]
subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(snps_sub))))
names(subhits)[6]<-"snps"
subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 

p1<-ggplot() + geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
                        stat="identity", color="gray37") +
  facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
  xlab("window") + ylab("sig. SNPs") +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank()) +
  xlab("position (Mbp)")
ggsave(paste0("~/Desktop/15_gwas_hits_", as.character(as.integer(win)), "DNA transposon_Unassigned.pdf"),
       p1, height=3, width=8)

### Retrotransposon Other
sub<-rpt[rpt$CLASS=="Retrotransposon" & rpt$FAMILY=="Other",]
snps_sub<-snps[snps$id %in% sub$ID,]
subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(snps_sub))))
names(subhits)[6]<-"snps"
subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 

p1<-ggplot() + geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
                        stat="identity", color="gray37") +
  facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
  xlab("window") + ylab("sig. SNPs") +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank()) +
  xlab("position (Mbp)")
ggsave(paste0("~/Desktop/15_gwas_hits_", as.character(as.integer(win)), "Retrotransposon_Unassigned.pdf"),
       p1, height=3, width=8)

