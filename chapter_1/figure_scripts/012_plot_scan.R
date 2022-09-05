# plot scan
# cjfiscus
# 2021-03-31

library(pacman)
p_load(ggplot2, cowplot, tidyr, matrixStats, ggpubr, Hmisc, reshape2, dplyr, 
       patchwork, ggrepel, viridis, patchwork)

# sR 100KB scan
df<-read.delim("~/Desktop/MANUSCRIPTS/athal_data/100KB.txt", check.names=F)

## transform
df[,2:ncol(df)]<-2^df[,2:ncol(df)]

## calculate sR
v<-as.data.frame(cbind(df$Feature, rowMedians(as.matrix(df[,2:ncol(df)])), 
                       rowMaxs(as.matrix(df[,2:ncol(df)])),
                       rowMins(as.matrix(df[,2:ncol(df)]))))

names(v)<-c("Feature", "median", "max", "min")

v$range<-as.numeric(as.character(v$max))-as.numeric(as.character(v$min))
v$rom<-as.numeric(as.character(v$range))/as.numeric(as.character(v$median))

v<-separate(v, Feature, c("chr","region"), sep=":")
v<-separate(v, region, c("start", "end"), sep="-")

v$start<-as.numeric(as.character(v$start))
v$end<-as.numeric(as.character(v$end))

v<-as.data.frame(v)
v$pos<-v$start + (v$end-v$start)/2

## add features
feats<-read.delim("~/Desktop/MANUSCRIPTS/athal_data/features.txt")

# standard scan
p1<-ggplot() + geom_bar(data=v, aes(x=pos, y=rom, group=chr), 
                    color="gray37", alpha=0.9, stat="identity") +
  facet_grid(.~chr, scales="free_x", space="free_x") +
  xlab("position (Mbp)") + ylab("sR") +
  theme_cowplot() +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank(),
    strip.text.x=element_text(face="bold", size=14),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.minor=element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  coord_cartesian(ylim=c(0,2.5)) +
  scale_y_continuous(expand = c(0, 0))
ggsave("scan.jpeg", p1, height=3, width=8)

##### scan with cents highlighted

scan<-ggplot() + 
  geom_rect(data=feats[feats$feat=="CEN",], 
            aes(xmin=start, xmax=end, ymin=0, ymax=5), fill="lightpink", 
            alpha=0.5) +
  geom_bar(data=v, aes(x=pos, y=rom, group=chr), 
                    color="gray37", alpha=0.9, stat="identity") +
  facet_grid(.~chr, scales="free_x", space="free_x") +
  xlab("position (Mbp)") + ylab("sR") +
  theme_cowplot() +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank(),
    strip.text.x=element_text(face="bold", size=14),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.minor=element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  coord_cartesian(ylim=c(0,2.5)) +
  scale_y_continuous(expand = c(0, 0))
ggsave("scan1.jpeg", scan, height=2, width=8)

#####
# generate list of top windows
h<-quantile(v$rom)[4] + 1.5*IQR(quantile(v$rom))
high<-v[v$rom > 0.18,]

# determine highly variable windows without cytological feature overlap
m<-merge(high, feats, by="chr")
m$overlap<-ifelse(m$start.x <= m$start.y & m$end.x >= m$start.y |
                    m$start.x >= m$start.y & m$end.x <= m$end.y | 
                    m$start.x <= m$end.y & m$end.x >= m$end.y,
                  "yes", "no")
m<-m[m$feat=="CEN",]
m<-as.data.frame(unique(m[,c("chr", "start.x", "end.x", "pos", "rom", "overlap")]))
m<-m[m$overlap=="no",]
write.table(m, "topwins.txt", sep="\t", quote=F, row.names=F)

p1<-ggplot() + 
  geom_rect(data=feats[feats$feat=="CEN",], 
            aes(xmin=start, xmax=end, ymin=0, ymax=5), fill="lightpink", 
            alpha=0.5) +
  geom_bar(data=v, aes(x=pos, y=rom, group=chr), 
           color="gray37", alpha=0.9, stat="identity") +
  geom_bar(data=m, aes(x=pos, y=rom, group=chr), 
           color="blue", alpha=0.9, stat="identity") +
  geom_hline(yintercept=0.18) +
  facet_grid(.~chr, scales="free_x", space="free_x") +
  xlab("position (Mbp)") + ylab("sR") +
  theme_cowplot() +
  scale_x_continuous(labels=function(x)x/1000000) +
  theme(
    strip.background = element_blank(),
    strip.text.x=element_text(face="bold", size=14),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.minor=element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  coord_cartesian(ylim=c(0,2.5)) +
  scale_y_continuous(expand = c(0, 0))
ggsave("scan_annot.jpeg", p1, height=3, width=8)

# distribution
p1<-ggplot(v, aes(x=rom)) + 
  geom_density() +
  geom_vline(aes(xintercept=h), color="red") +
  theme_cowplot()
ggsave("dist.jpeg", p1, height=3, width=3.5)

######
# gene and repeat scans
gene<-read.table("~/Desktop/MANUSCRIPTS/athal_data/cov_genes.txt")
names(gene)<-c("chr", "start", "end", "N", "bp", "total", "prop_genes")
gene$pos<-gene$start + (gene$end-gene$start)/2
gene<-gene[,c("chr","start","end", "pos", "prop_genes")]

repeats<-read.table("~/Desktop/MANUSCRIPTS/athal_data/cov_repeats.txt")
names(repeats)<-c("chr", "start", "end", "N", "bp", "total", "prop_repeats")
repeats$pos<-repeats$start + (repeats$end-repeats$start)/2
repeats<-repeats[,c("chr","start","end","pos", "prop_repeats")]

m<-merge(gene, v, by=c("chr", "pos"))
m<-merge(repeats, m, by=c("chr", "pos"))

p1<-ggplot(m, aes(x=prop_genes, y=log2(rom))) + 
  geom_point() + geom_smooth() +
  theme_cowplot() + 
  xlab("Proportion of genic sequence in bin") +
  ylab("log2(sR)")
ggsave("genic.jpeg", p1, height=3, width=4)

p1<-ggplot(m, aes(x=prop_repeats, y=log2(rom))) + 
  geom_point() + geom_smooth() +
  theme_cowplot() +
  xlab("Proportion of repetitive sequence in bin") +
  ylab("log2(sR)")
ggsave("rep.jpeg", p1, height=3, width=4)

## feature density strip plots
strip<-merge(gene, repeats, by=c("chr", "start", "end"))
#strip$prop_repeats<-scale(strip$prop_repeats, scale=F)
#strip$prop_genes<-scale(strip$prop_genes, scale=F)
strip<-strip[,c("chr", "start", "end", "prop_genes", "prop_repeats")]
names(strip)[4:5]<-c("gene", "repeat")
strip<-melt(strip, id.vars=c("chr", "start", "end"))
strip$pos<-(strip$start+strip$end)/2 + strip$start
strip$variable<-factor(strip$variable, levels=c("repeat", "gene"))

genescan<-ggplot(strip[strip$variable=="gene",], aes(x=pos, y=variable, fill=value, color=value)) +
  geom_raster(alpha=0.9, interpolate=T) +
  facet_grid(.~chr, scales="free_x", space="free_x") +
  theme_cowplot() +
  scale_fill_viridis(option="inferno", limits=c(0,1)) +
  scale_color_viridis(option="inferno", limits=c(0,1)) +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks=element_blank(),
    panel.grid.minor=element_blank(),
    axis.line=element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

repeatscan<-ggplot(strip[strip$variable=="repeat",], aes(x=pos, y=variable, fill=value, color=value)) +
  geom_raster(alpha=0.9, interpolate=T) +
  facet_grid(.~chr, scales="free_x", space="free_x") +
  theme_cowplot() +
  scale_fill_viridis(option="inferno", limits=c(0,1)) +
  scale_color_viridis(option="inferno", limits=c(0,1)) +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks=element_blank(),
    panel.grid.minor=element_blank(),
    axis.line=element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )
panel<-(genescan + theme(legend.position="none"))/
  (repeatscan + theme(legend.position="none"))
ggsave("den.jpeg", height=0.5, width=8)   

panel2<-scan / panel + plot_layout(heights=c(1, 0.2))
ggsave("all.pdf", height=3.5, width=12)
