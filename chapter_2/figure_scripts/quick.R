setwd("~/Desktop")
library(pacman)
p_load(ggplot2, cowplot, ggcorrplot, colorspace, tidyverse)

df<-read.delim("per_var_exp.txt")

p1<-ggplot(df, aes(x=pc, y=prop_var_exp)) + geom_bar(stat="identity") +
  theme_cowplot() +
  scale_x_continuous(breaks=seq(1,19,1)) + xlab("PC") + ylab("Percent variance explained")
ggsave("~/Desktop/pervar.pdf", p1, height=3, width=5)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/SCRIPTS")

df<-read.delim("../data/rad_data_samples_cw.txt")
names(df)[1]<-"sample"
df1<-read.delim("../data/bioclim_pc_scores.txt")
df<-merge(df, df1, by="sample")
df1<-read.delim("../data/cluster_assignments.txt")
df1<-df1[df1$K==2,]
df1<-as.data.frame(unique(df1[,c("sample", "assignment")]))
df<-merge(df, df1, by="sample")

lst<-c(13, 16, 19, 22, 23, 24)

cs<-data.frame()
for (i in lst){
  nm<-names(df)[i]
  c<-cor.test(df[,60], as.numeric(df[,i]), use="pairwise.complete.obs")
  c2<-cor.test(df[,61], as.numeric(df[,i]), use="pairwise.complete.obs")
  
  cs<-as.data.frame(rbind(cs, cbind("pc1", nm, c$p.value, c$estimate),
                          cbind("pc2", nm, c2$p.value, c2$estimate)))
}
cs$V4<-as.numeric(cs$V4)
names(cs)<-c("PC", "phenotype", "p", "R")

ord<-cs[cs$PC=="pc1",]
ord<-ord[order(ord$R),]
cs$phenotype<-factor(cs$phenotype, levels=ord$phenotype)

p1<-ggplot(cs, aes(x=PC, y=phenotype, fill=R)) + geom_tile() +
  scale_fill_gradient2(limits=c(-1,1)) + theme_cowplot()
ggsave("~/Desktop/cor.pdf", p1, height=4, width=6)


setwd("~/Desktop")
df<-read.delim("2022-07-25_CBP_GWAS_REGIONS - Sheet1.tsv")
df<-df[df$ID=="bioclim_PC1" | df$ID=="bioclim_PC2",]
df<-as.data.frame(unique(df[,c("ID", "CHR", "START", "END", "SIZE", "TAG", 
                               "MINOR", "MAJOR","MAF",
          "CANDIDATE", "NOTE")]))

df$id<-gsub("bioclim_", "", df$ID)
df$ID<-NULL

df1<-df %>% group_by(TAG) %>% summarize(ids=paste0(id, collapse=","), across()) 
df1$id<-NULL
df1<-as.data.frame(unique(df1))
df1$chr<-gsub("SCF_","",df1$CHR)
df1<-df1[order(as.numeric(df1$chr), as.numeric(df1$START)),]
df2<-df1[,c("ids", "TAG", "MINOR", "MAJOR", "MAF", "CHR", "START", "END")]
write.table(df2, "~/Desktop/table.txt", sep="\t", quote=F, row.names=F)

df3<-df1[nchar(df1$NOTE) > 0,]
df3<-df3[,c("TAG", "CANDIDATE", "NOTE")]
write.table(df3, "~/Desktop/table2.txt", sep="\t", quote=F, row.names=F)
