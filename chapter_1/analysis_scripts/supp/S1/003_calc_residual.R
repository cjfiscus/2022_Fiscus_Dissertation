#!/usr/bin/env Rscript

# load pkgs
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, ggplot2, matrixStats)

# load in data
df<-fread("../results/12mers.txt")

# linear model
model<-lm(log(`6909`+1) ~ log(TAIR10+1), data=df)
model

# write out residuals
df$model_residuals<-residuals(model)
write.Out<-as.data.frame(cbind(df$mer, df$model_residuals))
names(write.Out)<-c("mer", "residual")
write.table(write.Out, "../results/residuals.txt", sep="\t", quote=F, row.names=F)
