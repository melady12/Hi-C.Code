library(reshape2)
library(scales)
args<-commandArgs(TRUE)
gene<-args[1]
gefile<-args[2]
difile<-args[3]
tadfile<-args[4]
mat<-args[5]
pre<-args[6]
oup<-args[7]



gf<-read.table(gefile,stringsAsFactors=F,sep="\t")
gf<-gf[gf$V4==gene,]

med<-mean(c(gf$V2,gf$V3))
st<-floor((med-1500000)/25000)*25000
ed<-floor((med+1500000)/25000)*25000

lst<-st/25000+1
led<-ed/25000+1


MAT<-read.table(gzfile(mat),stringsAsFactors=F,sep="\t")
mdat<-MAT[lst:led,lst:led]
sname<-names(mdat)
mdat<-as.matrix(mdat)
mdat[lower.tri(mdat)]<-NA
colnames(mdat)<-1:length(sname)
rownames(mdat)<-1:length(sname)
#数据形式转化
pdat<-melt(mdat)
#去除NA值
pdat<-pdat[!is.na(pdat$value),]

#处理坐标
pdat$x<-(pdat$Var2+pdat$Var1)/sqrt(2)
pdat$y<-(pdat$Var2-pdat$Var1)/sqrt(2)

#pdat<-pdat[pdat$y<=c(max(pdat$y)/1.5),]


wid<-1*sqrt(2)/2
pdat$x1<-pdat$x-wid
pdat$x2<-pdat$x
pdat$x3<-pdat$x+wid
pdat$x4<-pdat$x
pdat$y1<-pdat$y
pdat$y2<-pdat$y-wid
pdat$y3<-pdat$y
pdat$y4<-pdat$y+wid
#加入NA，用于隔离点群
pdat$int<-NA
#进行颜色赋值
#quan<-quantile(pdat$value,probs=0.95)
quan<-100
pdat$value[pdat$value>=quan]<-quan
pdat$cols<-col_numeric(c('white','red'),domain=c(0,quan),na.color='gray')(pdat$value)
#polygon输入向量赋值
pdat<-pdat[pdat$y<=c(max(pdat$y)/1.5),]
whole.x<-as.vector(t(pdat[,c(6:9,14)]))
whole.y<-as.vector(t(pdat[,c(10:13,14)]))
#画板
pdf(oup,w=5,h=3)
layout(matrix(1:3,byrow=TRUE,nc=1),w=c(5,0),h=c(2.3,0.4,1.6))
#作图-高级绘图函数
par(xaxs='i',yaxs='i',tcl=-0.25,mgp=c(1.5,0.5,0),mar=c(1,5,1,2))
plot(1:10,type='n',axes=F,xlab='',ylab='',xlim=range(whole.x,na.rm=TRUE),ylim=range(whole.y,na.rm=TRUE)+c(wid,0),xaxs='i',yaxs='i',main=pre)
#作图-低级绘图函数polygon
polygon(x=whole.x,y=whole.y,col=pdat$cols,border=NA)
#作图-低级绘图函数axis-添加x轴
axis(1,at=range(whole.x,na.rm=TRUE),labels=paste(round(c(st,ed)/1000000,2),"Mb"),xpd=NA)
text(min(whole.x,na.rm=T),min(whole.y,na.rm=T),labels=gf$V1,xpd=NA,pos=2)
print(gf$V1)
gg<-rescale(c(st,gf$V2,gf$V3,ed),to=range(whole.x,na.rm=TRUE))
abline(v=gg[2:3],lty=2,lwd=0.5)
#legend

#par(fig=c(0.1,0.3,0.1,0.8),new=TRUE,mar=c(0,1,1,0))
#ran<-c(0,quan)
#ayy<-seq(ran[1],ran[2],len=100)
#cols<-col_numeric(c('white','red'),domain=c(0,quan))(ayy)
#lcols<-col_quantile('YlOrRd',ayy,n=12)(ayy)
#image(x=1:2,y=ayy,z=t(ayy),axes=F,xlab='',ylab='',xpd=NA,col=brewer.pal(9,'YlOrRd')[1:9])
#image(x=1:2,y=ayy,z=t(ayy),axes=F,xlab='',ylab='',xpd=NA,col=cols)
#mtext(2,text="interaction contacts)",line=0.5)
#mtext(2,text=expression(Log[2]~"("~Normalized),line=1.5)
#text(min(whole.x),0,pos=2,labels=gf$V1,xpd=NA)

#axis(4,las=2,cex.axis=0.9,lwd=0,lwd.ticks=1)
#box(bty='n')

#22222
tad<-read.table(tadfile,stringsAsFactors=F,sep="\t")
tad<-tad[tad$V1==gf$V1&tad$V2>=st&tad$V3<=ed,]
print(tad)
par(mar=c(0.1,5,1,2))
plot(1:10,type='n',axes=F,xlab='',ylab='',xlim=c(st,ed),ylim=c(0,1),xaxs='i',yaxs='i')
rect(tad$V2,0,tad$V3,1,col=c('#E1A949','#92B1BE'),border=NA)
text(st,0.5,pos=2,labels="TAD",xpd=NA)
abline(v=c(gf$V2,gf$V3),lty=2,lwd=0.5,xpd=NA)

#33333
di<-read.table(difile,stringsAsFactors=F,sep="\t")
di<-di[di$V1==gf$V1&di$V2>=st&di$V3<=ed,]
print(di)
di$cols<-ifelse(di$V4>0,"#a9252d","#2969a6")
par(mar=c(1,5,1,2))
plot(1:10,type='n',axes=F,xlab='',ylab='DI',xlim=c(st,ed),ylim=c(-60,60),xaxs='i',yaxs='i')
rect(di$V2,rep(0,length(di$V2)),di$V3,di$V4,col=di$cols,border=NA)
axis(2,las=2)
abline(v=c(gf$V2,gf$V3),lty=2,lwd=0.5,xpd=NA)

#legend
par(fig=c(0.1,0.15,0.6,0.9),new=TRUE,mar=c(0,1,1,0))
ran<-c(0,quan)
ayy<-seq(ran[1],ran[2],len=100)
cols<-col_numeric(c('white','red'),domain=c(0,quan))(ayy)
#lcols<-col_quantile('YlOrRd',ayy,n=12)(ayy)
#image(x=1:2,y=ayy,z=t(ayy),axes=F,xlab='',ylab='',xpd=NA,col=brewer.pal(9,'YlOrRd')[1:9])
image(x=1:2,y=ayy,z=t(ayy),axes=F,xlab='',ylab='',xpd=NA,col=cols)
mtext(2,text="interaction contacts)",line=0.5,cex=0.7)
mtext(2,text=expression(Log[2]~"("~Normalized),line=1.5,cex=0.7)
text(min(whole.x),0,pos=2,labels=gf$V1,xpd=NA)

axis(4,las=2,cex.axis=0.9,lwd=0,lwd.ticks=1,cex.axis=0.8)
box(bty='n')


dev.off()


