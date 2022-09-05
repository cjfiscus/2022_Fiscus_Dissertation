#!/usr/bin/env Rscript
# correct p vals in prep for meta
# cjfiscus
# 2021-07-16

args = commandArgs(trailingOnly=TRUE)

library(data.table)

## import p vals and chr
filename=args[1]

df<-fread(filename)

## calculate lambda
df$chisq<-qchisq(df$p_lrt, 1)
l<-median(df$chisq)/qchisq(0.5,1)

# correct p values on non-focal chrs, focal chr p = 1
df$p_corrected<-df$p_lrt/l

fwrite(df, paste0(filename,"2"), sep="\t", quote=F, row.names=F)
