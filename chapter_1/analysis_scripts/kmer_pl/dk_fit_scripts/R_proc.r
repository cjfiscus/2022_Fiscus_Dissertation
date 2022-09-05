#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
##########################################################################
#libs
library(pacman)
#
adjthese<-function(x,y){
  themod<-lm(x~y)
  return(themod$residuals+themod$coefficients[1])
}
#
p_load(data.table, ggplot2, cowplot, rsvd)
##########################################################################
if (length(args)<4) {
  stop("Expecting four arguements acc.file samplehead kmer.mat out.file.n", call.=FALSE)
}
##########################################################################
#Read supplemental data
acc<-read.delim(args[1])
rownames(acc)<-acc[,1]
thehead<-read.delim(args[2],header=F,stringsAsFactors=F)
##########################################################################
df<-as.data.frame(fread(args[3],stringsAsFactors = F,header=F))
colnames(df)<-thehead
dfmer<-df$mer
df$mer<-NULL
##########################################################################
#hardcode outgroups
outgroups<-c("A_lyrata1",
             "A_lyrata2",
             "C_rubella1",
             "C_rubella2")
acc<-acc[colnames(df[,!colnames(df)%in%outgroups]),]
df.o<-df[,!colnames(df)%in%rownames(acc)]
df.at<-df[,colnames(df)%in%rownames(acc)]

#
corrected<-apply(as.matrix(df.at),1,adjthese,acc$seq_by)
out<-cbind(dfmer,t(corrected),df.o)
write.table(out,args[4],quote=F,row.names=F,col.names=F,sep="\t")
