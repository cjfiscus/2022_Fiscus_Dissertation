#!/usr/bin/env Rscript

# install/ load required libs 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, matrixStats, stringr)

# correct raw kmer counts by proportion of lib
# in each gc bin 

# read sysargs
args = commandArgs(trailingOnly=TRUE)

## FUNCTIONS ## 
# define GC correction function 
CorrectGC<-function(column){
  ## aggregate bins 
  gc<-data.table(cbind(GC_content, column))
  names(gc)<-c("GC_content", "column")
  gc<-as.data.frame(gc[,sum(column),by=GC_content])
  names(gc)<-c("GC_content", "column")
  
  ## merge reference table with table for column
  gc<-merge(REF, gc, by="GC_content")
  
  ## calculate weights
  gc$Proportion2<-gc$column/as.numeric(sum(gc$column))
  weights<-gc$Proportion/gc$Proportion2
  names(weights)<-gc$GC_content
  
  # column name that is being processed
  gc$Library<-names(counts)[i]
  i<<-i+1 # this is bad practice but could not find another way
  
  proportions<<-rbind(proportions, cbind(gc$Library, gc$GC_content, round(gc$Proportion2, 4)))
  
  ## apply weights
  Correction<-as.data.frame(cbind(column, GC_content))
  keys<-sapply(Correction[,2], function(x) weights[[as.character(x)]]) # get factors 
  Correction$Cts<-round(Correction$column*keys) # apply correction
  return(Correction$Cts)
}

#####

# read in K-mer counts 
counts<-fread(args[1])
dim(counts)
counts[1:5,1:5]

# Calculate row medians  
choose_cols=colnames(counts)[2:length(counts)] # skip first col
counts[, med:=rowMedians(as.matrix(.SD)),.SDcols=choose_cols]

# calculate GC content
GC_content<-str_count(counts$mer, pattern="G|C")

# aggregate countsums by gc 
meds<-as.data.table(cbind(GC_content, counts$med))
names(meds)<-c("GC_content", "Count")

REF<-meds[,sum(Count, na.rm=T), by=GC_content]
names(REF)<-c("GC_content", "Sum")
REF$Proportion<-REF$Sum/as.numeric(sum(REF$Sum))

write.table(REF, args[3], sep="\t", row.names=T, quote=F)

# cleanup 
counts$med<-NULL
rm(meds)

# data table to collect proportions from each library (calculated during correction)
proportions<-data.frame(Library=character(), GC_content=character(), Proportion=character())

# correct counts 
print("correction")
i<-2
counts[, (choose_cols):=lapply(.SD, CorrectGC), .SDcols=choose_cols]

# write kmer table out
fwrite(counts, file=args[2], sep="\t", quote=F)
rm(counts)

# write proportions table out
names(proportions)<-c("Library", "GC_content", "Proportion")
fwrite(proportions, file=args[4], sep="\t", quote=F)

sessionInfo()
