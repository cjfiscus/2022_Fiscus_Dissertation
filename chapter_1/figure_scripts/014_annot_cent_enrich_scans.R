#!/usr/bin/env Rscript
# windowed SNPs
# cjfiscus
# 2021-12-16

library(pacman)
p_load(data.table, IRanges, GenomicRanges, ggplot2, cowplot)

## define cent enriched or not
df1<-read.delim("~/Desktop/data/fishers_incent_gwas.txt")
enrich<-df1[df1$p < 0.05/nrow(df1),]
noenrich<-df1[!df1$ID %in% enrich$ID,]

## repeat annots
annot<-read.table("~/Desktop/data/TAIR10_repeats_custom.bed")
names(annot)<-c("CHR", "START", "END", "ID")
annot$LEN<-annot$END-annot$START

### filter by id
annot<-annot[annot$ID %in% df1$ID,]

### filter by len
annot<-annot[annot$LEN >= 50,]

### filter orgs
annot<-annot[!annot$CHR=="Mt",]
annot<-annot[!annot$CHR=="Pt",]
#########

## save annot to GRanges
rngs<-GRanges(seqnames = annot$CHR, 
              ranges=IRanges(start=as.numeric(annot$START), end=as.numeric(annot$END)), id=annot$ID)

## set seq lens
ranges<-read.delim("~/Desktop/data/TAIR10_ranges.txt")
ranges<-ranges[1:5,]
rng<-ranges$END
names(rng)<-ranges$CHR
seqlengths(rngs)<-rng

## count hits per window
### define window
win<-100000
bins<-tileGenome(seqlengths(rngs), tilewidth=win, cut.last.tile.in.chrom=TRUE)

### Cent enriched
rngs_sub<-rngs[rngs$id %in% enrich$ID,]
subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(rngs_sub))))
names(subhits)[6]<-"snps"
subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 

feats<-read.delim("~/Desktop/data/features.txt")
names(feats)[1]<-"seqnames"

p1<-ggplot() +
  geom_rect(data=feats[feats$feat=="CEN",], 
            aes(xmin=start, xmax=end, ymin=0, ymax=Inf), fill="lightpink", 
            alpha=0.5) +
  geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
                        stat="identity", color="gray37") +
  facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
  xlab("window") + ylab("annotations per bin") +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank()) +
  xlab("position (Mbp)")
ggsave(paste0("~/Desktop/14_cent_enrich_annot_", as.character(as.integer(win)), ".pdf"),
       p1, height=3, width=8)

### Not cent enriched
rngs_sub<-rngs[rngs$id %in% noenrich$ID,]
subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(rngs_sub))))
names(subhits)[6]<-"snps"
subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 

feats<-read.delim("~/Desktop/data/features.txt")
names(feats)[1]<-"seqnames"

p1<-ggplot() +
  geom_rect(data=feats[feats$feat=="CEN",], 
            aes(xmin=start, xmax=end, ymin=0, ymax=Inf), fill="lightpink", 
            alpha=0.5) +
  geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
           stat="identity", color="gray37") +
  facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
  xlab("window") + ylab("annotations per bin") +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank()) +
  xlab("position (Mbp)")
ggsave(paste0("~/Desktop/14_cent_notenrich_annot_", as.character(as.integer(win)), ".pdf"),
       p1, height=3, width=8)


### All considered
rngs_sub<-rngs[rngs$id %in% df1$ID,]
subhits<-as.data.frame(cbind(as.data.frame(bins), countOverlaps(bins, unique(rngs_sub))))
names(subhits)[6]<-"snps"
subhits$pos<-subhits$start+(subhits$end-subhits$start)/2 

feats<-read.delim("~/Desktop/data/features.txt")
names(feats)[1]<-"seqnames"

p1<-ggplot() +
  geom_rect(data=feats[feats$feat=="CEN",], 
            aes(xmin=start, xmax=end, ymin=0, ymax=Inf), fill="lightpink", 
            alpha=0.5) +
  geom_bar(data=subhits, aes(x=pos, y=snps, group=seqnames), 
           stat="identity", color="gray37") +
  facet_grid(.~seqnames, scales="free_x", space="free_x") + theme_cowplot() +
  xlab("window") + ylab("annotations per bin") +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank()) +
  xlab("position (Mbp)")
ggsave(paste0("~/Desktop/14_cent_alltested_annot_", as.character(as.integer(win)), ".pdf"),
       p1, height=3, width=8)
