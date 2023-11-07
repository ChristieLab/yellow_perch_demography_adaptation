#========================================================================================#
# Script created by Xiaoshen Yin, all rights reserved, contact at xiaoshenyin@hotmail.com
# Script created in version R 4.2.0
# This script: calculate observed heterozygosity
#========================================================================================#

## set directory
setwd("/scratch/bell/yin168/pflavescens/gt_w210_01/fig1d_het")

## calculate observed heterozygosity with window size of 5 MB and step size of 2.5 MB
dat<-read.table("ypsnp_greenbay_mindp03_ind090_bi_maf005.GT.FORMAT",head=T)
dat$g00<-rowSums(dat[,3:dim(dat)[2]]=="0|0"|dat[,3:dim(dat)[2]]=="0/0")
dat$g01<-rowSums(dat[,3:dim(dat)[2]]=="0|1"|dat[,3:dim(dat)[2]]=="0/1")
dat$g10<-rowSums(dat[,3:dim(dat)[2]]=="1|0"|dat[,3:dim(dat)[2]]=="1/0")
dat$g11<-rowSums(dat[,3:dim(dat)[2]]=="1|1"|dat[,3:dim(dat)[2]]=="1/1")
dat$missing<-rowSums(dat[,3:dim(dat)[2]]=="."|dat[,3:dim(dat)[2]]=="./.")
dat$total<-(dat$g00+dat$g01+dat$g10+dat$g11+dat$missing)
unique(dat$total)
dat$homo<-(dat$g00+dat$g11)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$het<-(dat$g01+dat$g10)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$totalfreq<-dat$homo+dat$het
unique(dat$totalfreq)
write.table(dat,"het_greenbay_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
chr<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
posmean<-NULL
for (i in 1:length(chr)){
posmeani<-NULL
for (j in 1:(max(dat[dat$CHROM==chr[i],]$POS)/2500000)){
meanj<-mean(dat[(dat$CHROM==chr[i]&dat$POS>=((j-1)*5000000-(j-1)*2500000+1)&dat$POS<=((j*5000000)-(j-1)*2500000)),]$het)
posmeanj<-data.frame(chr[i],((j-1)*5000000-(j-1)*2500000+1),meanj)
posmeani<-rbind(posmeani,posmeanj)
}
posmean<-rbind(posmean,posmeani)
}
write.table(posmean,"meanhet_win5000kstep2500k_greenbay_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")

dat<-read.table("ypsnp_muskegon_mindp03_ind090_bi_maf005.GT.FORMAT",head=T)
dat$g00<-rowSums(dat[,3:dim(dat)[2]]=="0|0"|dat[,3:dim(dat)[2]]=="0/0")
dat$g01<-rowSums(dat[,3:dim(dat)[2]]=="0|1"|dat[,3:dim(dat)[2]]=="0/1")
dat$g10<-rowSums(dat[,3:dim(dat)[2]]=="1|0"|dat[,3:dim(dat)[2]]=="1/0")
dat$g11<-rowSums(dat[,3:dim(dat)[2]]=="1|1"|dat[,3:dim(dat)[2]]=="1/1")
dat$missing<-rowSums(dat[,3:dim(dat)[2]]=="."|dat[,3:dim(dat)[2]]=="./.")
dat$total<-(dat$g00+dat$g01+dat$g10+dat$g11+dat$missing)
unique(dat$total)
dat$homo<-(dat$g00+dat$g11)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$het<-(dat$g01+dat$g10)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$totalfreq<-dat$homo+dat$het
unique(dat$totalfreq)
write.table(dat,"het_muskegon_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
chr<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
posmean<-NULL
for (i in 1:length(chr)){
posmeani<-NULL
for (j in 1:(max(dat[dat$CHROM==chr[i],]$POS)/2500000)){
meanj<-mean(dat[(dat$CHROM==chr[i]&dat$POS>=((j-1)*5000000-(j-1)*2500000+1)&dat$POS<=((j*5000000)-(j-1)*2500000)),]$het)
posmeanj<-data.frame(chr[i],((j-1)*5000000-(j-1)*2500000+1),meanj)
posmeani<-rbind(posmeani,posmeanj)
}
posmean<-rbind(posmean,posmeani)
}
write.table(posmean,"meanhet_win5000kstep2500k_muskegon_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")

dat<-read.table("ypsnp_grandhaven_mindp03_ind090_bi_maf005.GT.FORMAT",head=T)
dat$g00<-rowSums(dat[,3:dim(dat)[2]]=="0|0"|dat[,3:dim(dat)[2]]=="0/0")
dat$g01<-rowSums(dat[,3:dim(dat)[2]]=="0|1"|dat[,3:dim(dat)[2]]=="0/1")
dat$g10<-rowSums(dat[,3:dim(dat)[2]]=="1|0"|dat[,3:dim(dat)[2]]=="1/0")
dat$g11<-rowSums(dat[,3:dim(dat)[2]]=="1|1"|dat[,3:dim(dat)[2]]=="1/1")
dat$missing<-rowSums(dat[,3:dim(dat)[2]]=="."|dat[,3:dim(dat)[2]]=="./.")
dat$total<-(dat$g00+dat$g01+dat$g10+dat$g11+dat$missing)
unique(dat$total)
dat$homo<-(dat$g00+dat$g11)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$het<-(dat$g01+dat$g10)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$totalfreq<-dat$homo+dat$het
unique(dat$totalfreq)
write.table(dat,"het_grandhaven_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
chr<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
posmean<-NULL
for (i in 1:length(chr)){
posmeani<-NULL
for (j in 1:(max(dat[dat$CHROM==chr[i],]$POS)/2500000)){
meanj<-mean(dat[(dat$CHROM==chr[i]&dat$POS>=((j-1)*5000000-(j-1)*2500000+1)&dat$POS<=((j*5000000)-(j-1)*2500000)),]$het)
posmeanj<-data.frame(chr[i],((j-1)*5000000-(j-1)*2500000+1),meanj)
posmeani<-rbind(posmeani,posmeanj)
}
posmean<-rbind(posmean,posmeani)
}
write.table(posmean,"meanhet_win5000kstep2500k_grandhaven_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")

dat<-read.table("ypsnp_michigancity_mindp03_ind090_bi_maf005.GT.FORMAT",head=T)
dat$g00<-rowSums(dat[,3:dim(dat)[2]]=="0|0"|dat[,3:dim(dat)[2]]=="0/0")
dat$g01<-rowSums(dat[,3:dim(dat)[2]]=="0|1"|dat[,3:dim(dat)[2]]=="0/1")
dat$g10<-rowSums(dat[,3:dim(dat)[2]]=="1|0"|dat[,3:dim(dat)[2]]=="1/0")
dat$g11<-rowSums(dat[,3:dim(dat)[2]]=="1|1"|dat[,3:dim(dat)[2]]=="1/1")
dat$missing<-rowSums(dat[,3:dim(dat)[2]]=="."|dat[,3:dim(dat)[2]]=="./.")
dat$total<-(dat$g00+dat$g01+dat$g10+dat$g11+dat$missing)
unique(dat$total)
dat$homo<-(dat$g00+dat$g11)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$het<-(dat$g01+dat$g10)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$totalfreq<-dat$homo+dat$het
unique(dat$totalfreq)
write.table(dat,"het_michigancity_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
chr<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
posmean<-NULL
for (i in 1:length(chr)){
posmeani<-NULL
for (j in 1:(max(dat[dat$CHROM==chr[i],]$POS)/2500000)){
meanj<-mean(dat[(dat$CHROM==chr[i]&dat$POS>=((j-1)*5000000-(j-1)*2500000+1)&dat$POS<=((j*5000000)-(j-1)*2500000)),]$het)
posmeanj<-data.frame(chr[i],((j-1)*5000000-(j-1)*2500000+1),meanj)
posmeani<-rbind(posmeani,posmeanj)
}
posmean<-rbind(posmean,posmeani)
}
write.table(posmean,"meanhet_win5000kstep2500k_michigancity_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")

dat<-read.table("ypsnp_milwaukee_mindp03_ind090_bi_maf005.GT.FORMAT",head=T)
dat$g00<-rowSums(dat[,3:dim(dat)[2]]=="0|0"|dat[,3:dim(dat)[2]]=="0/0")
dat$g01<-rowSums(dat[,3:dim(dat)[2]]=="0|1"|dat[,3:dim(dat)[2]]=="0/1")
dat$g10<-rowSums(dat[,3:dim(dat)[2]]=="1|0"|dat[,3:dim(dat)[2]]=="1/0")
dat$g11<-rowSums(dat[,3:dim(dat)[2]]=="1|1"|dat[,3:dim(dat)[2]]=="1/1")
dat$missing<-rowSums(dat[,3:dim(dat)[2]]=="."|dat[,3:dim(dat)[2]]=="./.")
dat$total<-(dat$g00+dat$g01+dat$g10+dat$g11+dat$missing)
unique(dat$total)
dat$homo<-(dat$g00+dat$g11)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$het<-(dat$g01+dat$g10)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$totalfreq<-dat$homo+dat$het
unique(dat$totalfreq)
write.table(dat,"het_milwaukee_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
chr<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
posmean<-NULL
for (i in 1:length(chr)){
posmeani<-NULL
for (j in 1:(max(dat[dat$CHROM==chr[i],]$POS)/2500000)){
meanj<-mean(dat[(dat$CHROM==chr[i]&dat$POS>=((j-1)*5000000-(j-1)*2500000+1)&dat$POS<=((j*5000000)-(j-1)*2500000)),]$het)
posmeanj<-data.frame(chr[i],((j-1)*5000000-(j-1)*2500000+1),meanj)
posmeani<-rbind(posmeani,posmeanj)
}
posmean<-rbind(posmean,posmeani)
}
write.table(posmean,"meanhet_win5000kstep2500k_milwaukee_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")

dat<-read.table("ypsnp_naubinway_mindp03_ind090_bi_maf005.GT.FORMAT",head=T)
dat$g00<-rowSums(dat[,3:dim(dat)[2]]=="0|0"|dat[,3:dim(dat)[2]]=="0/0")
dat$g01<-rowSums(dat[,3:dim(dat)[2]]=="0|1"|dat[,3:dim(dat)[2]]=="0/1")
dat$g10<-rowSums(dat[,3:dim(dat)[2]]=="1|0"|dat[,3:dim(dat)[2]]=="1/0")
dat$g11<-rowSums(dat[,3:dim(dat)[2]]=="1|1"|dat[,3:dim(dat)[2]]=="1/1")
dat$missing<-rowSums(dat[,3:dim(dat)[2]]=="."|dat[,3:dim(dat)[2]]=="./.")
dat$total<-(dat$g00+dat$g01+dat$g10+dat$g11+dat$missing)
unique(dat$total)
dat$homo<-(dat$g00+dat$g11)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$het<-(dat$g01+dat$g10)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$totalfreq<-dat$homo+dat$het
unique(dat$totalfreq)
write.table(dat,"het_naubinway_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
chr<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
posmean<-NULL
for (i in 1:length(chr)){
posmeani<-NULL
for (j in 1:(max(dat[dat$CHROM==chr[i],]$POS)/2500000)){
meanj<-mean(dat[(dat$CHROM==chr[i]&dat$POS>=((j-1)*5000000-(j-1)*2500000+1)&dat$POS<=((j*5000000)-(j-1)*2500000)),]$het)
posmeanj<-data.frame(chr[i],((j-1)*5000000-(j-1)*2500000+1),meanj)
posmeani<-rbind(posmeani,posmeanj)
}
posmean<-rbind(posmean,posmeani)
}
write.table(posmean,"meanhet_win5000kstep2500k_naubinway_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")

dat<-read.table("ypsnp_suttonsbay_mindp03_ind090_bi_maf005.GT.FORMAT",head=T)
dat$g00<-rowSums(dat[,3:dim(dat)[2]]=="0|0"|dat[,3:dim(dat)[2]]=="0/0")
dat$g01<-rowSums(dat[,3:dim(dat)[2]]=="0|1"|dat[,3:dim(dat)[2]]=="0/1")
dat$g10<-rowSums(dat[,3:dim(dat)[2]]=="1|0"|dat[,3:dim(dat)[2]]=="1/0")
dat$g11<-rowSums(dat[,3:dim(dat)[2]]=="1|1"|dat[,3:dim(dat)[2]]=="1/1")
dat$missing<-rowSums(dat[,3:dim(dat)[2]]=="."|dat[,3:dim(dat)[2]]=="./.")
dat$total<-(dat$g00+dat$g01+dat$g10+dat$g11+dat$missing)
unique(dat$total)
dat$homo<-(dat$g00+dat$g11)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$het<-(dat$g01+dat$g10)/(dat$g00+dat$g01+dat$g10+dat$g11)
dat$totalfreq<-dat$homo+dat$het
unique(dat$totalfreq)
write.table(dat,"het_suttonsbay_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
chr<-c("CM014935.1","CM014936.1","CM014937.1","CM014938.1","CM014939.1","CM014940.1","CM014941.1","CM014942.1","CM014943.1","CM014944.1","CM014945.1","CM014946.1","CM014947.1","CM014948.1","CM014949.1","CM014950.1","CM014951.1","CM014952.1","CM014953.1","CM014954.1","CM014955.1","CM014956.1","CM014957.1","CM014958.1")
posmean<-NULL
for (i in 1:length(chr)){
posmeani<-NULL
for (j in 1:(max(dat[dat$CHROM==chr[i],]$POS)/2500000)){
meanj<-mean(dat[(dat$CHROM==chr[i]&dat$POS>=((j-1)*5000000-(j-1)*2500000+1)&dat$POS<=((j*5000000)-(j-1)*2500000)),]$het)
posmeanj<-data.frame(chr[i],((j-1)*5000000-(j-1)*2500000+1),meanj)
posmeani<-rbind(posmeani,posmeanj)
}
posmean<-rbind(posmean,posmeani)
}
write.table(posmean,"meanhet_win5000kstep2500k_suttonsbay_mindp03_ind090_bi_maf005.txt",row.names=F,quote=F,sep="\t")
