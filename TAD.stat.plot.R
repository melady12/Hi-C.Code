library(ggplot2)
library(plyr)
library(gtools)
library(RColorBrewer)
args<-commandArgs(TRUE)
outdir<-args[1]
infile<-paste(args[1],"/total.combined.domaincalls",sep="")
outfile.table<-paste(args[1],"/TAD.bychr.csv",sep="")
outfile.pdf<-paste(args[1],"/TAD.plot.pdf",sep="")
outfile.png<-paste(args[1],"/TAD.plot.png",sep="")
data<-read.table(infile,stringsAsFactors=F,sep="\t")
names(data)=c('Chr','Start','End')
inter2<-data
inter2$Chr<-"Total"
#inter2<-rbind(data,inter2)
out1<-ddply(data,.(Chr),summarize,Number=length(Chr),Length.average=mean(End-Start+1),Length.median=median(End-Start+1))
out2<-ddply(inter2,.(Chr),summarize,Number=length(Chr),Length.average=mean(End-Start+1),Length.median=median(End-Start+1))
out1<-out1[mixedorder(out1$Chr),]
outres<-rbind(out1,out2)
write.table(outres,outfile.table,quote=F,sep=",",row.names=F)

inter<-data
inter$Chr<-"zzgenome"
dat<-rbind(data,inter)
dat$len<-(dat$End-dat$Start+1)/1000000
dat<-dat[mixedorder(dat$Chr),]
row.names(dat)<-1:dim(dat)[1]
dat$Chr<-factor(dat$Chr,levels=unique(dat$Chr))
cols<-rep("white",length(unique(dat$Chr)))
cols[length(cols)]<-brewer.pal(9,"Set1")[2]
label<-as.character(unique(dat$Chr))
label[label=="zzgenome"]<-"genome"
breaks=seq(floor(min(dat$len)),ceiling(max(dat$len)),2)
p<-ggplot(dat,aes(Chr,len))+geom_boxplot(fill=cols,outlier.size = 1.5,)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_text(angl=45,hjust=1))+scale_y_continuous(breaks=breaks)+scale_x_discrete(labels=label)+xlab("Chromosome")+ylab("Size in Mb")

figwid<-10*length(unique(dat$Chr))/24
fig.height<-6.04*ceiling(max(dat$len))/14

ggsave(outfile.png, plot=p,type="cairo-png", width=figwid, height=fig.height, dpi=700)
ggsave(outfile.pdf, plot=p, width=figwid, height=fig.height)

