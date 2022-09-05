#!/usr/bin/env Rscript
# prep eGWAS
# cjfiscus
# 2021-12-07
# edited 2022-07-09

library(pacman)
p_load(readxl, raster, rgdal)

# import data and slice out relevant info
df<-read_excel("../data/IITA-Cowpea collection.xls", sheet="Accession passport data")
df<-as.data.frame(df)
df<-df[,c("Accession name", "Latitude", "Longitude", "Collecting/acquisition source", "Biological status of accession")]

## subset to accessions we have
acc<-read.table("../data/samples.txt")
df<-df[df$`Accession name` %in% acc$V1,]

## remove market collected
rm<-"Market or shop"
df1<-df[!df$`Collecting/acquisition source` %in% rm,]
rm<-"Breeding / Research material"
df1<-df1[!df1$`Biological status of accession` %in% rm,]
df1<-df1[,1:3]
df1<-na.omit(df1)
df1<-df1[,c("Accession name", "Longitude", "Latitude")]

# compile wc data
files<-list.files(path="../data", pattern=".tif", full.names=T)
rasterStack<-stack(files)

row.names(df1)<-df1$`Accession name`
df1$`Accession name`<-NULL

worldclimData<-as.data.frame(cbind(row.names(df1), extract(rasterStack, df1)))
names(worldclimData)<-c("ID", "BIO1", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", 
                        "BIO15", "BIO16", "BIO17", "BIO18", "BIO19", "BIO2", 
                        "BIO3", "BIO4", "BIO5", "BIO6", "BIO7", "BIO8", "BIO9")
worldclimData<-worldclimData[,c("ID", paste0("BIO", seq(1,19)))]

## merge with coordinates
names(df)[1]<-"ID"
m<-merge(df, worldclimData, by="ID")

write.table(worldclimData, "../data/cowpea_worldclim2_1.txt", sep="\t", quote=F, row.names=F)
#####

# write out data as tfam 
names(acc)<-"ID"

m<-merge(acc, worldclimData, by="ID", all.x=T)
m<-m[match(acc$ID, m$ID),]

info<-as.data.frame(cbind(m$ID, m$ID, 0, 0, 0))
m<-as.data.frame(cbind(info, m[,2:ncol(m)]))

write.table(m, "../data/cowpea_envgwas.fam", sep=" ", quote=F, row.names=F, col.names=F)

## export climate vars
vars<-names(m)[6:ncol(m)]
write.table(vars, "../data/worldclim_vars.txt", sep="\t", quote=F, row.names=F, col.names=F)
