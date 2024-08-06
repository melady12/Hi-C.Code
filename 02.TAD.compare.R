args<-commandArgs(TRUE)
library(TADCompare)
library(dplyr)
library(SpectralTAD)
library(GenomicRanges)


BH<-args[1]
M2<-args[2]
chr<-args[3]
oup<-args[4]



bh<-read.table(gzfile(BH),stringsAsFactors=F,sep="\t")
bh<-as.matrix(bh)

m2<-read.table(gzfile(M2),stringsAsFactors=F,sep="\t")
m2<-as.matrix(m2)


tbh<-read.table("./BH.TAD.bed",stringsAsFactors=F,sep="\t")
tm2<-read.table("./M2.TAD.bed",stringsAsFactors=F,sep="\t")

cbh<-tbh[tbh$V1==chr,]
cm2<-tm2[tm2$V1==chr,]
names(cbh)<-c('chr','start','end')
names(cm2)<-c('chr','start','end')
Combined_Bed = list(as_tibble(cbh), as_tibble(cm2))
#Combined_Bed = list(GRanges(cbh), GRanges(cm2))
colnames(bh)<-(1:(dim(bh)[2])-1)*20000
rownames(bh)<-1:(dim(bh)[2])
colnames(m2)<-(1:(dim(m2)[2])-1)*20000
rownames(m2)<-1:(dim(m2)[2])
TD_Compare <-  TADCompare(bh, m2, resolution = 20000, pre_tads = Combined_Bed)$TAD_Frame[,1:5]
TD_Compare$chr<-chr
TD_Compare<-TD_Compare[,c(6,1:5)]

write.table(TD_Compare,oup,quote=F,sep="\t",row.names=F)
