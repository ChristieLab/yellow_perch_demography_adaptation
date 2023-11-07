#========================================================================================#
# Script created by Xiaoshen Yin, all rights reserved, contact at xiaoshenyin@hotmail.com
# Script created in version R 4.2.0
# This script: make Tajima's D plots
#========================================================================================#

## set directory
setwd("E:/Purdue/perch/MS/figure/figure2/fig2c")

## calculate mean Tajima's D with window size of 5 MB
muskegon_dat1<-read.table("muskegon_maf002_5000k.Tajima.D",head=T)
muskegon_dat1[,3]<-NULL
colnames(muskegon_dat1)<-c("CHROM","POS","TajimaD")
muskegon_dat2<-muskegon_dat1
muskegon_dat3<-muskegon_dat2[muskegon_dat2$CHROM=="CM014935.1"|muskegon_dat2$CHROM=="CM014936.1"|muskegon_dat2$CHROM=="CM014937.1"|muskegon_dat2$CHROM=="CM014938.1"|muskegon_dat2$CHROM=="CM014939.1"|muskegon_dat2$CHROM=="CM014940.1"|muskegon_dat2$CHROM=="CM014941.1"|muskegon_dat2$CHROM=="CM014942.1"|muskegon_dat2$CHROM=="CM014943.1"|muskegon_dat2$CHROM=="CM014944.1"|muskegon_dat2$CHROM=="CM014945.1"|muskegon_dat2$CHROM=="CM014946.1"|muskegon_dat2$CHROM=="CM014947.1"|muskegon_dat2$CHROM=="CM014948.1"|muskegon_dat2$CHROM=="CM014949.1"|muskegon_dat2$CHROM=="CM014950.1"|muskegon_dat2$CHROM=="CM014951.1"|muskegon_dat2$CHROM=="CM014952.1"|muskegon_dat2$CHROM=="CM014953.1"|muskegon_dat2$CHROM=="CM014954.1"|muskegon_dat2$CHROM=="CM014955.1"|muskegon_dat2$CHROM=="CM014956.1"|muskegon_dat2$CHROM=="CM014957.1"|muskegon_dat2$CHROM=="CM014958.1",]
chr1<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1")
max.pos<-sapply(chr1,function(i){(max(muskegon_dat3[muskegon_dat3$CHROM==i,]$POS))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
added.length<-(add1+add2)
names(added.length) <-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
muskegon_dat3$gp<-muskegon_dat3$POS+added.length[muskegon_dat3$CHROM]
muskegon_dat3$color<-"#FFC000"

greenbay_dat1<-read.table("greenbay_maf002_5000k.Tajima.D",head=T)
greenbay_dat1[,3]<-NULL
colnames(greenbay_dat1)<-c("CHROM","POS","TajimaD")
greenbay_dat2<-greenbay_dat1
greenbay_dat3<-greenbay_dat2[greenbay_dat2$CHROM=="CM014935.1"|greenbay_dat2$CHROM=="CM014936.1"|greenbay_dat2$CHROM=="CM014937.1"|greenbay_dat2$CHROM=="CM014938.1"|greenbay_dat2$CHROM=="CM014939.1"|greenbay_dat2$CHROM=="CM014940.1"|greenbay_dat2$CHROM=="CM014941.1"|greenbay_dat2$CHROM=="CM014942.1"|greenbay_dat2$CHROM=="CM014943.1"|greenbay_dat2$CHROM=="CM014944.1"|greenbay_dat2$CHROM=="CM014945.1"|greenbay_dat2$CHROM=="CM014946.1"|greenbay_dat2$CHROM=="CM014947.1"|greenbay_dat2$CHROM=="CM014948.1"|greenbay_dat2$CHROM=="CM014949.1"|greenbay_dat2$CHROM=="CM014950.1"|greenbay_dat2$CHROM=="CM014951.1"|greenbay_dat2$CHROM=="CM014952.1"|greenbay_dat2$CHROM=="CM014953.1"|greenbay_dat2$CHROM=="CM014954.1"|greenbay_dat2$CHROM=="CM014955.1"|greenbay_dat2$CHROM=="CM014956.1"|greenbay_dat2$CHROM=="CM014957.1"|greenbay_dat2$CHROM=="CM014958.1",]
chr1<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1")
max.pos<-sapply(chr1,function(i){(max(greenbay_dat3[greenbay_dat3$CHROM==i,]$POS))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
added.length<-(add1+add2)
names(added.length) <-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
greenbay_dat3$gp<-greenbay_dat3$POS+added.length[greenbay_dat3$CHROM]
greenbay_dat3$color<-"#FF0000"

grandhaven_dat1<-read.table("grandhaven_maf002_5000k.Tajima.D",head=T)
grandhaven_dat1[,3]<-NULL
colnames(grandhaven_dat1)<-c("CHROM","POS","TajimaD")
grandhaven_dat2<-grandhaven_dat1
grandhaven_dat3<-grandhaven_dat2[grandhaven_dat2$CHROM=="CM014935.1"|grandhaven_dat2$CHROM=="CM014936.1"|grandhaven_dat2$CHROM=="CM014937.1"|grandhaven_dat2$CHROM=="CM014938.1"|grandhaven_dat2$CHROM=="CM014939.1"|grandhaven_dat2$CHROM=="CM014940.1"|grandhaven_dat2$CHROM=="CM014941.1"|grandhaven_dat2$CHROM=="CM014942.1"|grandhaven_dat2$CHROM=="CM014943.1"|grandhaven_dat2$CHROM=="CM014944.1"|grandhaven_dat2$CHROM=="CM014945.1"|grandhaven_dat2$CHROM=="CM014946.1"|grandhaven_dat2$CHROM=="CM014947.1"|grandhaven_dat2$CHROM=="CM014948.1"|grandhaven_dat2$CHROM=="CM014949.1"|grandhaven_dat2$CHROM=="CM014950.1"|grandhaven_dat2$CHROM=="CM014951.1"|grandhaven_dat2$CHROM=="CM014952.1"|grandhaven_dat2$CHROM=="CM014953.1"|grandhaven_dat2$CHROM=="CM014954.1"|grandhaven_dat2$CHROM=="CM014955.1"|grandhaven_dat2$CHROM=="CM014956.1"|grandhaven_dat2$CHROM=="CM014957.1"|grandhaven_dat2$CHROM=="CM014958.1",]
chr1<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1")
max.pos<-sapply(chr1,function(i){(max(grandhaven_dat3[grandhaven_dat3$CHROM==i,]$POS))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
added.length<-(add1+add2)
names(added.length) <-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
grandhaven_dat3$gp<-grandhaven_dat3$POS+added.length[grandhaven_dat3$CHROM]
grandhaven_dat3$color<-"#92D050"

michigancity_dat1<-read.table("michigancity_maf002_5000k.Tajima.D",head=T)
michigancity_dat1[,3]<-NULL
colnames(michigancity_dat1)<-c("CHROM","POS","TajimaD")
michigancity_dat2<-michigancity_dat1
michigancity_dat3<-michigancity_dat2[michigancity_dat2$CHROM=="CM014935.1"|michigancity_dat2$CHROM=="CM014936.1"|michigancity_dat2$CHROM=="CM014937.1"|michigancity_dat2$CHROM=="CM014938.1"|michigancity_dat2$CHROM=="CM014939.1"|michigancity_dat2$CHROM=="CM014940.1"|michigancity_dat2$CHROM=="CM014941.1"|michigancity_dat2$CHROM=="CM014942.1"|michigancity_dat2$CHROM=="CM014943.1"|michigancity_dat2$CHROM=="CM014944.1"|michigancity_dat2$CHROM=="CM014945.1"|michigancity_dat2$CHROM=="CM014946.1"|michigancity_dat2$CHROM=="CM014947.1"|michigancity_dat2$CHROM=="CM014948.1"|michigancity_dat2$CHROM=="CM014949.1"|michigancity_dat2$CHROM=="CM014950.1"|michigancity_dat2$CHROM=="CM014951.1"|michigancity_dat2$CHROM=="CM014952.1"|michigancity_dat2$CHROM=="CM014953.1"|michigancity_dat2$CHROM=="CM014954.1"|michigancity_dat2$CHROM=="CM014955.1"|michigancity_dat2$CHROM=="CM014956.1"|michigancity_dat2$CHROM=="CM014957.1"|michigancity_dat2$CHROM=="CM014958.1",]
chr1<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1")
max.pos<-sapply(chr1,function(i){(max(michigancity_dat3[michigancity_dat3$CHROM==i,]$POS))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
added.length<-(add1+add2)
names(added.length) <-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
michigancity_dat3$gp<-michigancity_dat3$POS+added.length[michigancity_dat3$CHROM]
michigancity_dat3$color<-"#00B0F0"

milwaukee_dat1<-read.table("milwaukee_maf002_5000k.Tajima.D",head=T)
milwaukee_dat1[,3]<-NULL
colnames(milwaukee_dat1)<-c("CHROM","POS","TajimaD")
milwaukee_dat2<-milwaukee_dat1
milwaukee_dat3<-milwaukee_dat2[milwaukee_dat2$CHROM=="CM014935.1"|milwaukee_dat2$CHROM=="CM014936.1"|milwaukee_dat2$CHROM=="CM014937.1"|milwaukee_dat2$CHROM=="CM014938.1"|milwaukee_dat2$CHROM=="CM014939.1"|milwaukee_dat2$CHROM=="CM014940.1"|milwaukee_dat2$CHROM=="CM014941.1"|milwaukee_dat2$CHROM=="CM014942.1"|milwaukee_dat2$CHROM=="CM014943.1"|milwaukee_dat2$CHROM=="CM014944.1"|milwaukee_dat2$CHROM=="CM014945.1"|milwaukee_dat2$CHROM=="CM014946.1"|milwaukee_dat2$CHROM=="CM014947.1"|milwaukee_dat2$CHROM=="CM014948.1"|milwaukee_dat2$CHROM=="CM014949.1"|milwaukee_dat2$CHROM=="CM014950.1"|milwaukee_dat2$CHROM=="CM014951.1"|milwaukee_dat2$CHROM=="CM014952.1"|milwaukee_dat2$CHROM=="CM014953.1"|milwaukee_dat2$CHROM=="CM014954.1"|milwaukee_dat2$CHROM=="CM014955.1"|milwaukee_dat2$CHROM=="CM014956.1"|milwaukee_dat2$CHROM=="CM014957.1"|milwaukee_dat2$CHROM=="CM014958.1",]
chr1<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1")
max.pos<-sapply(chr1,function(i){(max(milwaukee_dat3[milwaukee_dat3$CHROM==i,]$POS))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
added.length<-(add1+add2)
names(added.length) <-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
milwaukee_dat3$gp<-milwaukee_dat3$POS+added.length[milwaukee_dat3$CHROM]
milwaukee_dat3$color<-"#FF66FF"

suttonsbay_dat1<-read.table("suttonsbay_maf002_5000k.Tajima.D",head=T)
suttonsbay_dat1[,3]<-NULL
colnames(suttonsbay_dat1)<-c("CHROM","POS","TajimaD")
suttonsbay_dat2<-suttonsbay_dat1
suttonsbay_dat3<-suttonsbay_dat2[suttonsbay_dat2$CHROM=="CM014935.1"|suttonsbay_dat2$CHROM=="CM014936.1"|suttonsbay_dat2$CHROM=="CM014937.1"|suttonsbay_dat2$CHROM=="CM014938.1"|suttonsbay_dat2$CHROM=="CM014939.1"|suttonsbay_dat2$CHROM=="CM014940.1"|suttonsbay_dat2$CHROM=="CM014941.1"|suttonsbay_dat2$CHROM=="CM014942.1"|suttonsbay_dat2$CHROM=="CM014943.1"|suttonsbay_dat2$CHROM=="CM014944.1"|suttonsbay_dat2$CHROM=="CM014945.1"|suttonsbay_dat2$CHROM=="CM014946.1"|suttonsbay_dat2$CHROM=="CM014947.1"|suttonsbay_dat2$CHROM=="CM014948.1"|suttonsbay_dat2$CHROM=="CM014949.1"|suttonsbay_dat2$CHROM=="CM014950.1"|suttonsbay_dat2$CHROM=="CM014951.1"|suttonsbay_dat2$CHROM=="CM014952.1"|suttonsbay_dat2$CHROM=="CM014953.1"|suttonsbay_dat2$CHROM=="CM014954.1"|suttonsbay_dat2$CHROM=="CM014955.1"|suttonsbay_dat2$CHROM=="CM014956.1"|suttonsbay_dat2$CHROM=="CM014957.1"|suttonsbay_dat2$CHROM=="CM014958.1",]
chr1<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1")
max.pos<-sapply(chr1,function(i){(max(suttonsbay_dat3[suttonsbay_dat3$CHROM==i,]$POS))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
added.length<-(add1+add2)
names(added.length) <-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
suttonsbay_dat3$gp<-suttonsbay_dat3$POS+added.length[suttonsbay_dat3$CHROM]
suttonsbay_dat3$color<-"#7030A0"

naubinway_dat1<-read.table("naubinway_maf002_5000k.Tajima.D",head=T)
naubinway_dat1[,3]<-NULL
colnames(naubinway_dat1)<-c("CHROM","POS","TajimaD")
naubinway_dat2<-naubinway_dat1
naubinway_dat3<-naubinway_dat2[naubinway_dat2$CHROM=="CM014935.1"|naubinway_dat2$CHROM=="CM014936.1"|naubinway_dat2$CHROM=="CM014937.1"|naubinway_dat2$CHROM=="CM014938.1"|naubinway_dat2$CHROM=="CM014939.1"|naubinway_dat2$CHROM=="CM014940.1"|naubinway_dat2$CHROM=="CM014941.1"|naubinway_dat2$CHROM=="CM014942.1"|naubinway_dat2$CHROM=="CM014943.1"|naubinway_dat2$CHROM=="CM014944.1"|naubinway_dat2$CHROM=="CM014945.1"|naubinway_dat2$CHROM=="CM014946.1"|naubinway_dat2$CHROM=="CM014947.1"|naubinway_dat2$CHROM=="CM014948.1"|naubinway_dat2$CHROM=="CM014949.1"|naubinway_dat2$CHROM=="CM014950.1"|naubinway_dat2$CHROM=="CM014951.1"|naubinway_dat2$CHROM=="CM014952.1"|naubinway_dat2$CHROM=="CM014953.1"|naubinway_dat2$CHROM=="CM014954.1"|naubinway_dat2$CHROM=="CM014955.1"|naubinway_dat2$CHROM=="CM014956.1"|naubinway_dat2$CHROM=="CM014957.1"|naubinway_dat2$CHROM=="CM014958.1",]
chr1<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1")
max.pos<-sapply(chr1,function(i){(max(naubinway_dat3[naubinway_dat3$CHROM==i,]$POS))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
added.length<-(add1+add2)
names(added.length) <-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
naubinway_dat3$gp<-naubinway_dat3$POS+added.length[naubinway_dat3$CHROM]
naubinway_dat3$color<-"#0000CC"

## make Tajima's D plots
pdf("TajimaD_by_win_5000k_maf002.pdf",width=9,height=3)
par(mar=c(1,1,1,1))
plot(muskegon_dat3[muskegon_dat3$CHROM=="CM014935.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014935.1",]$TajimaD,xlab="Genomic position (bp)",ylab="Tajima's D",xlim=c(0,795000000),ylim=c(-2,2),col=muskegon_dat3[muskegon_dat3$CHROM=="CM014935.1",]$color,type="l",axes=F,xaxs="i",yaxs="i")
axis(side=1,at=c(0,1000000000),lwd=1.5,lwd.ticks=0)
# axis(side=1,at=c(0,200000000,400000000,600000000,800000000),lwd=1.5,lwd.ticks=1.5)
axis(side=2,at=c(-2,-1,0,1,2),lwd=1.5,lwd.ticks=1.5)
# abline(h=0.1,lwd=1.5)
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014936.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014936.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014936.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014937.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014937.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014937.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014938.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014938.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014938.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014939.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014939.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014939.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014940.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014940.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014940.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014941.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014941.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014941.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014942.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014942.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014942.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014943.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014943.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014943.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014944.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014944.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014944.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014945.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014945.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014945.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014946.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014946.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014946.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014947.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014947.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014947.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014948.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014948.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014948.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014949.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014949.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014949.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014950.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014950.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014950.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014951.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014951.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014951.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014952.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014952.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014952.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014953.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014953.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014953.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014954.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014954.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014954.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014955.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014955.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014955.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014956.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014956.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014956.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014957.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014957.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014957.1",]$color,type="l")
points(muskegon_dat3[muskegon_dat3$CHROM=="CM014958.1",]$gp,muskegon_dat3[muskegon_dat3$CHROM=="CM014958.1",]$TajimaD,col=muskegon_dat3[muskegon_dat3$CHROM=="CM014958.1",]$color,type="l")

points(greenbay_dat3[greenbay_dat3$CHROM=="CM014935.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014935.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014935.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014936.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014936.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014936.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014937.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014937.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014937.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014938.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014938.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014938.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014939.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014939.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014939.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014940.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014940.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014940.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014941.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014941.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014941.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014942.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014942.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014942.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014943.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014943.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014943.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014944.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014944.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014944.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014945.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014945.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014945.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014946.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014946.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014946.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014947.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014947.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014947.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014948.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014948.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014948.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014949.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014949.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014949.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014950.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014950.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014950.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014951.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014951.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014951.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014952.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014952.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014952.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014953.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014953.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014953.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014954.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014954.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014954.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014955.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014955.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014955.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014956.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014956.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014956.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014957.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014957.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014957.1",]$color,type="l")
points(greenbay_dat3[greenbay_dat3$CHROM=="CM014958.1",]$gp,greenbay_dat3[greenbay_dat3$CHROM=="CM014958.1",]$TajimaD,col=greenbay_dat3[greenbay_dat3$CHROM=="CM014958.1",]$color,type="l")

points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014935.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014935.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014935.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014936.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014936.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014936.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014937.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014937.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014937.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014938.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014938.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014938.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014939.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014939.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014939.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014940.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014940.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014940.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014941.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014941.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014941.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014942.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014942.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014942.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014943.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014943.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014943.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014944.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014944.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014944.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014945.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014945.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014945.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014946.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014946.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014946.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014947.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014947.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014947.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014948.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014948.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014948.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014949.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014949.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014949.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014950.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014950.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014950.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014951.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014951.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014951.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014952.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014952.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014952.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014953.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014953.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014953.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014954.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014954.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014954.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014955.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014955.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014955.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014956.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014956.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014956.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014957.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014957.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014957.1",]$color,type="l")
points(grandhaven_dat3[grandhaven_dat3$CHROM=="CM014958.1",]$gp,grandhaven_dat3[grandhaven_dat3$CHROM=="CM014958.1",]$TajimaD,col=grandhaven_dat3[grandhaven_dat3$CHROM=="CM014958.1",]$color,type="l")

points(michigancity_dat3[michigancity_dat3$CHROM=="CM014935.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014935.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014935.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014936.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014936.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014936.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014937.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014937.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014937.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014938.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014938.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014938.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014939.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014939.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014939.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014940.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014940.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014940.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014941.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014941.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014941.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014942.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014942.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014942.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014943.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014943.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014943.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014944.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014944.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014944.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014945.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014945.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014945.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014946.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014946.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014946.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014947.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014947.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014947.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014948.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014948.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014948.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014949.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014949.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014949.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014950.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014950.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014950.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014951.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014951.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014951.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014952.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014952.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014952.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014953.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014953.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014953.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014954.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014954.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014954.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014955.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014955.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014955.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014956.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014956.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014956.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014957.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014957.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014957.1",]$color,type="l")
points(michigancity_dat3[michigancity_dat3$CHROM=="CM014958.1",]$gp,michigancity_dat3[michigancity_dat3$CHROM=="CM014958.1",]$TajimaD,col=michigancity_dat3[michigancity_dat3$CHROM=="CM014958.1",]$color,type="l")

points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014935.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014935.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014935.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014936.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014936.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014936.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014937.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014937.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014937.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014938.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014938.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014938.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014939.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014939.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014939.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014940.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014940.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014940.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014941.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014941.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014941.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014942.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014942.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014942.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014943.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014943.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014943.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014944.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014944.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014944.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014945.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014945.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014945.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014946.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014946.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014946.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014947.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014947.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014947.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014948.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014948.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014948.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014949.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014949.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014949.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014950.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014950.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014950.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014951.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014951.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014951.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014952.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014952.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014952.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014953.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014953.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014953.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014954.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014954.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014954.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014955.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014955.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014955.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014956.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014956.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014956.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014957.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014957.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014957.1",]$color,type="l")
points(milwaukee_dat3[milwaukee_dat3$CHROM=="CM014958.1",]$gp,milwaukee_dat3[milwaukee_dat3$CHROM=="CM014958.1",]$TajimaD,col=milwaukee_dat3[milwaukee_dat3$CHROM=="CM014958.1",]$color,type="l")

points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014935.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014935.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014935.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014936.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014936.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014936.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014937.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014937.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014937.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014938.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014938.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014938.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014939.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014939.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014939.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014940.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014940.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014940.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014941.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014941.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014941.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014942.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014942.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014942.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014943.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014943.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014943.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014944.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014944.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014944.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014945.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014945.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014945.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014946.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014946.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014946.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014947.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014947.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014947.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014948.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014948.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014948.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014949.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014949.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014949.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014950.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014950.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014950.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014951.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014951.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014951.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014952.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014952.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014952.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014953.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014953.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014953.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014954.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014954.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014954.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014955.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014955.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014955.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014956.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014956.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014956.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014957.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014957.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014957.1",]$color,type="l")
points(suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014958.1",]$gp,suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014958.1",]$TajimaD,col=suttonsbay_dat3[suttonsbay_dat3$CHROM=="CM014958.1",]$color,type="l")

points(naubinway_dat3[naubinway_dat3$CHROM=="CM014935.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014935.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014935.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014936.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014936.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014936.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014937.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014937.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014937.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014938.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014938.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014938.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014939.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014939.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014939.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014940.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014940.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014940.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014941.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014941.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014941.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014942.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014942.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014942.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014943.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014943.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014943.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014944.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014944.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014944.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014945.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014945.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014945.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014946.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014946.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014946.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014947.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014947.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014947.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014948.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014948.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014948.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014949.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014949.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014949.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014950.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014950.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014950.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014951.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014951.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014951.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014952.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014952.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014952.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014953.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014953.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014953.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014954.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014954.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014954.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014955.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014955.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014955.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014956.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014956.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014956.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014957.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014957.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014957.1",]$color,type="l")
points(naubinway_dat3[naubinway_dat3$CHROM=="CM014958.1",]$gp,naubinway_dat3[naubinway_dat3$CHROM=="CM014958.1",]$TajimaD,col=naubinway_dat3[naubinway_dat3$CHROM=="CM014958.1",]$color,type="l")

dev.off()
