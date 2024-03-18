library(scales)
dat<-read.table("genes.for.plot.txt",stringsAsFactors=F,sep="\t")
Mi<-read.table("BH.20Kb.PC1.xls",stringsAsFactors=F)
FMi<-read.table("M2.20Kb.PC1.xls",stringsAsFactors=F)

pdf("AB.swtich.new.pdf",w=4.5,h=2.5)
for (i in seq(dat$V1)){
#for (i in 1:2){
	med<-mean(c(dat$V2[i],dat$V3[i]))
	st<-floor((med-1500000)/25000)*25000
	ed<-ceiling((med+1500000)/25000)*25000
	pd1<-Mi[Mi$V1==dat$V1[i]&Mi$V2>=st&Mi$V3<=ed,]
	pd1$cols<-ifelse(pd1$V4>0,"#DA486D","#2385A9")
	pd2<-FMi[FMi$V1==dat$V1[i]&FMi$V2>=st&FMi$V3<=ed,]
	pd2$cols<-ifelse(pd2$V4>0,"#DA486D","#2385A9")
	ymt<-range(c(pd1$V4,pd2$V4))
	print(range(ymt))
	#layout(matrix(1:2,nc=1),w=c(5,0),h=c(1.5,1.5))
	par(fig=c(0,1,0,1),xaxs='i',yaxs='i',xpd=NA,lend=2,mar=c(3,2,2,1.5),mgp=c(1.3,0.4,0),new=FALSE,tcl=-0.15)
	plot(1:10,type='n',axes=F,xlab="Chromosomal position (Mb)",ylab='',xlim=c(st-800000,ed),main=substitute(italic(P),list(P=dat$V5[i])))
	box(bty='o')
	axis(1,at=c(st,ed),labels=round(c(st,ed)/1000000,2),xpd=NA)
	text(x=rep(st-800000,2),y=c(3,8.5),labels=c('BH','M2'),xpd=NA,pos=2)
	rect(floor(dat$V2[i]/25000)*25000,1,ceiling(dat$V3[i]/25000)*25000,10,border='black',lwd=0.5,lty=3,xpd=NA)
	par(fig=c(0,1,0.3,0.5),new=TRUE,mar=c(0,2,0,1.5))
	plot(1:10,type='n',axes=F,xlab='',ylab='',xlim=c(st-50000,ed),ylim=ymt)
	rect(pd1$V2,0,pd1$V3,pd1$V4,border=NA,col=pd1$cols)
	text(x=st+25000,y=0,labels="B",pos=1,xpd=NA)
	text(x=st+25000,y=0,labels="A",pos=3,xpd=NA)
	par(fig=c(0,1,0.55,0.75),new=TRUE,mar=c(0,2,0,1.5))
        plot(1:10,type='n',axes=F,xlab='',ylab='',xlim=c(st-50000,ed),ylim=ymt)
        rect(pd2$V2,0,pd2$V3,pd2$V4,border=NA,col=pd2$cols)
        text(x=st+25000,y=0,labels="B",pos=1,xpd=NA)
        text(x=st+25000,y=0,labels="A",pos=3,xpd=NA)
}
dev.off()
