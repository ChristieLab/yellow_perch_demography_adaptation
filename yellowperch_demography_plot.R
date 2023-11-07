#========================================================================================#
# Script created by Xiaoshen Yin, all rights reserved, contact at xiaoshenyin@hotmail.com
# Script created in version R 4.2.0
# This script: make demography plots
#========================================================================================#

## set directory
setwd("E:/Purdue/perch/MS/figure/figure3/fig3de/output")

## install packages
install.packages("scales",repos="https://mirror-hk.koddos.net/CRAN/")

library(scales)

## read in data and calculate median effective population size
greenbay<-NULL
for (j in 1:1000){
greenbayi<-read.table(paste0("select_ypsnp_greenbay_mindp03_ind090_bi_maf005_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt"))
greenbay<-data.frame(rbind(greenbay,greenbayi))
}

greenbaymedian<-NULL
for (i in 1:max(greenbay$V1)){
greenbaymediani<-cbind(i,median(greenbay[greenbay$V1==i,]$V2))
greenbaymedian<-data.frame(rbind(greenbaymedian,greenbaymediani))
}

muskegon<-NULL
for (j in 1:1000){
muskegoni<-read.table(paste0("select_ypsnp_muskegon_mindp03_ind090_bi_maf005_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt"))
muskegon<-data.frame(rbind(muskegon,muskegoni))
}

muskegonmedian<-NULL
for (i in 1:max(muskegon$V1)){
muskegonmediani<-cbind(i,median(muskegon[muskegon$V1==i,]$V2))
muskegonmedian<-data.frame(rbind(muskegonmedian,muskegonmediani))
}

xmax<-max(max(greenbaymedian$i),max(muskegonmedian$i))
ymax<-max(max(greenbaymedian$V2),max(muskegonmedian$V2))

## make demography plots
pdf("E:/Purdue/perch/MS/figure/figure3/fig3de/fig3d_y750_median_greenbay.pdf",height=5,width=5)
par(mar=c(1,1,1,1))
plot(0,0,xlim=c(0,xmax),ylim=c(0,200),col="#FFFFFF",axes=F,xaxs="i",yaxs="i")
for (j in 1:1000){
greenbayi<-read.table(paste0("select_ypsnp_greenbay_mindp03_ind090_bi_maf005_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt"))
lines(greenbayi$V1,greenbayi$V2,lwd=1.5,col=alpha("#FF0000",0.1))}
points(greenbaymedian$i,greenbaymedian$V2,type="l",lwd=1.5,col="#000000")
axis(side=1,at=c(0,150,300,450,600,750),lwd=1.5,lwd.ticks=1.5)
axis(side=2,at=c(0,50,100,150,200),lwd=1.5,lwd.ticks=1.5)
dev.off()

pdf("E:/Purdue/perch/MS/figure/figure3/fig3de/fig3d_y200_median_greenbay.pdf",height=5,width=5)
par(mar=c(1,1,1,1))
plot(0,0,xlim=c(0,200),ylim=c(0,200),col="#FFFFFF",axes=F,xaxs="i",yaxs="i")
for (j in 1:1000){
greenbayi<-read.table(paste0("select_ypsnp_greenbay_mindp03_ind090_bi_maf005_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt"))
lines(greenbayi$V1,greenbayi$V2,lwd=1.5,col=alpha("#FF0000",0.1))}
points(greenbaymedian$i,greenbaymedian$V2,type="l",lwd=1.5,col="#000000")
axis(side=1,at=c(0,50,100,150,200),lwd=1.5,lwd.ticks=1.5)
axis(side=2,at=c(0,50,100,150,200),lwd=1.5,lwd.ticks=1.5)
dev.off()

pdf("E:/Purdue/perch/MS/figure/figure3/fig3de/fig3e_y750_median_muskegon.pdf",height=5,width=5)
par(mar=c(1,1,1,1))
plot(0,0,xlim=c(0,xmax),ylim=c(0,200),col="#FFFFFF",axes=F,xaxs="i",yaxs="i")
for (j in 1:1000){
muskegoni<-read.table(paste0("select_ypsnp_muskegon_mindp03_ind090_bi_maf005_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt"))
lines(muskegoni$V1,muskegoni$V2,lwd=1.5,col=alpha("#FFC000",0.1))}
points(muskegonmedian$i,muskegonmedian$V2,type="l",lwd=1.5,col="#000000")
axis(side=1,at=c(0,150,300,450,600,750),lwd=1.5,lwd.ticks=1.5)
axis(side=2,at=c(0,50,100,150,200),lwd=1.5,lwd.ticks=1.5)
dev.off()

pdf("E:/Purdue/perch/MS/figure/figure3/fig3de/fig3e_y200_median_muskegon.pdf",height=5,width=5)
par(mar=c(1,1,1,1))
plot(0,0,xlim=c(0,200),ylim=c(0,200),col="#FFFFFF",axes=F,xaxs="i",yaxs="i")
for (j in 1:1000){
muskegoni<-read.table(paste0("select_ypsnp_muskegon_mindp03_ind090_bi_maf005_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt"))
lines(muskegoni$V1,muskegoni$V2,lwd=1.5,col=alpha("#FFC000",0.1))}
points(muskegonmedian$i,muskegonmedian$V2,type="l",lwd=1.5,col="#000000")
axis(side=1,at=c(0,50,100,150,200),lwd=1.5,lwd.ticks=1.5)
axis(side=2,at=c(0,50,100,150,200),lwd=1.5,lwd.ticks=1.5)
dev.off()
