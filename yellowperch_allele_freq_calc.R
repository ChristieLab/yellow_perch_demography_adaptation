#============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com
# Script created in version R 4.1.0 on 11/12/2021
# This script: takes vcf files (one_per_population) and calculates Ho and MAF quickly (~5 minutes or 2 million snps)
# Usage notes: if interested in differences in # snps among populations, would want to call snps within each population (e.g., 2x SNPs in Greenbay than Mainbaisin)
#============================================================================================================#
## Set working directory, import packages, source functions, initialize global variables
install.packages("gaston",repos="https://mirror-hk.koddos.net/CRAN/")
library(gaston)
setwd("/scratch/bell/yin168/pflavescens/gt_w210/figx_maf")

list.files()

d9<-read.vcf("ypsnp_mindp03_ind090_bi_muskegon.recode.vcf", convert.chr = FALSE)
gh18<-read.vcf("ypsnp_mindp03_ind090_bi_grandhaven.recode.vcf", convert.chr = FALSE)
gbtu18<-read.vcf("ypsnp_mindp03_ind090_bi_greenbay.recode.vcf", convert.chr = FALSE)
imc18<-read.vcf("ypsnp_mindp03_ind090_bi_michigancity.recode.vcf", convert.chr = FALSE)
mil19<-read.vcf("ypsnp_mindp03_ind090_bi_milwaukee.recode.vcf", convert.chr = FALSE)
mtc18<-read.vcf("ypsnp_mindp03_ind090_bi_suttonsbay.recode.vcf", convert.chr = FALSE)
nub18<-read.vcf("ypsnp_mindp03_ind090_bi_naubinway.recode.vcf", convert.chr = FALSE)

## set directory
setwd("/scratch/bell/yin168/pflavescens/gt_w210/figx_maf/allfreqdist")

## Start isolating populations
d9<-as.matrix(d9)
gh18<-as.matrix(gh18)
gbtu18<-as.matrix(gbtu18)
imc18<-as.matrix(imc18)
mil19<-as.matrix(mil19)
mtc18<-as.matrix(mtc18)
nub18<-as.matrix(nub18)

populations <- list()
populations[[1]] <- d9
populations[[2]] <- gh18
populations[[3]] <- gbtu18
populations[[4]] <- imc18
populations[[5]] <- mil19
populations[[6]] <- mtc18
populations[[7]] <- nub18

mafs <- list()
hos <- list()

for(x in 1:length(populations)){
  pop <- populations[[x]]
  
  ids  <- rownames(pop)
  print(ids) 
  pops <- substring(ids, 1, 3)
  table(pops)
  
  counts <- apply(pop, 2, table) # count across all rows
  sums   <- lapply(counts, sum)  # sum table counts

  MAJOR <- function(counts, n) {
    counts[names(counts)=="0"] * 2   # 2 because 2 alleles from homozygote, if statements for missing
    }
  
  HET <- function(counts, n) {
    counts[names(counts)=="1"]    # * nothing  because 1 allele from homozygote
    }
  
  # calculate total count of reference allele
  major   <- lapply(counts, MAJOR)  
  m1      <- which(lapply(major, length) == 0)  # add 0s if not found
  major[m1] <- 0
   
  # calculate total count of heterozygote reference allele
  het     <- lapply(counts, HET)
  m1      <- which(lapply(het, length) == 0) # add 0s if not found
  het[m1] <- 0
  
  # calculate minor allele frequency
  maf <- unlist(major) + unlist(het)
  maf <- maf/(unlist(sums) * 2)  # times 2 because 2 alleles per genotype
  maf <- 1-maf
  
  # some reference alleles may actually be the minor allele, if so calculate true minor allele freq
  m2 <- which(maf >= 0.5)
  maf[m2] <- 1-maf[m2]
   
  # name list elements for output 
  mafs[[x]] <- maf   
  names(mafs)[x] <- row.names(populations[[x]])[1]   # name elements in list
  
  # calculate observed heterozygosity
  ho <- unlist(het)/unlist(sums)
  hos[[x]] <- ho   
  names(hos)[x] <- row.names(populations[[x]])[1]   # name elements in list
  
}

## write.table(OUT, "mafs_and_Ho.txt", col.names = TRUE, sep="\t", append = FALSE)
  
for(n in 1:length(mafs)){
  
  pop <- mafs[[n]]
  #print(names(mafs[1]))
  #print(length(pop))
  hist(pop, breaks = 30, xlim = c(-0.05, 0.55), main = paste(names(mafs[n]), length(pop), sep = "_"), xlab = "minor allele frequence")
  abline(h=0)
  abline(v=0, col="blue")
  out1 <- mean(pop)
  out2 <- median(pop)
  abline(v=out1, col="red")
  abline(v=out2, col="orange")
  print(out1)
  print(out2)

}

for(n in 1:length(mafs)){
  
  pop <- hos[[n]]
  #print(names(mafs[1]))
  #print(length(pop))
  hist(pop, breaks = 30, xlim = c(-0.05, 1.05), main = paste(names(mafs[n]), length(pop), sep = "_"), xlab = "observed heterozygosity")
  abline(h=0)
  abline(v=0, col="blue")
  out1 <- mean(pop)
  out2 <- median(pop)
  abline(v=out1, col="red")
  abline(v=out2, col="orange")
  print(out1)
  print(out2)
  
}


temp <- data.frame(unlist(mafs))
nms  <- names(unlist(mafs))

populations[[1]] <- d9
populations[[2]] <- gh18
populations[[3]] <- gbtu18
populations[[4]] <- imc18
populations[[5]] <- mil19
populations[[6]] <- mtc18
populations[[7]] <- nub18

dat1 <- cbind(data.frame("d9", mafs[[1]]))
dat2 <- cbind(data.frame("gh", mafs[[2]]))
dat3 <- cbind(data.frame("gbt", mafs[[3]]))
dat4 <- cbind(data.frame("imc", mafs[[4]]))
dat5 <- cbind(data.frame("mil", mafs[[5]]))
dat6 <- cbind(data.frame("mtc", mafs[[6]]))
dat7 <- cbind(data.frame("nub", mafs[[7]]))

colnames(dat1) <- c("population", "mafs")
colnames(dat2) <- c("population", "mafs")
colnames(dat3) <- c("population", "mafs")
colnames(dat4) <- c("population", "mafs")
colnames(dat5) <- c("population", "mafs")
colnames(dat6) <- c("population", "mafs")
colnames(dat7) <- c("population", "mafs")

output <- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7)
head(output)
tail(output)
write.table(output, "mafs_maf000.txt", col.names = TRUE, sep="\t", append = FALSE)
