#=====================================================================================#
# Script created by Peter Euclide, all rights reserved, contact at peuclide@purdue.edu
# This script: takes MAF files to make MAF density plots
#=====================================================================================#

## set working directory, import packages, source functions, initialize global variables
install.packages("ggplot2")
install.packages("ggridges")

library(ggplot2)
library(ggridges)

## set RGB code
Muskegon     <-  rgb(255, 192, 0, max = 255)
GreenBay     <-  rgb(255, 0, 0, max = 255)
GrandHaven   <-  rgb(146, 208, 80, max = 255)
MichiganCity <-  rgb(0, 176, 240, max = 255)
SuttonsBay   <-  rgb(112, 48, 160, max = 255)
Naubinway    <-  rgb(0, 0, 204, max = 255)
Milwaukee    <-  rgb(255, 102, 255, max = 255)

## set directory
setwd("E:/Purdue/perch/figure/figure2/MAFs")

## load data
dat <- read.table("MAFs.txt")

## make MAF density plots
dat$population <- as.factor(dat$population)
levels(dat$population) <- c("Muskegon", "Green Bay", "Grand Haven", "Sutton's Bay", "Milwaukee", "Michigan City", "Naubinway")
dat$population <- factor(dat$population, c("Green Bay","Muskegon",  "Grand Haven", "Michigan City", "Milwaukee",  "Naubinway","Sutton's Bay"))
colnames(dat) <- c("Population", "Minor allele frequency")

p6<-ggplot(dat, aes(x = MAF, y = Population, fill = Population, color = Population)) +
  geom_density_ridges(scale = 5, bandwidth = .015, alpha = 0.5,quantile_lines = TRUE, quantiles = 2,size=0.6)+
  scale_y_discrete(limits=rev)+
  scale_x_continuous("Minor Allele Frequency", breaks = c(.1, .2, .3, .4, .5))+
  scale_fill_manual(values = c(GreenBay, Muskegon, GrandHaven, MichiganCity, Milwaukee, Naubinway, SuttonsBay))+
  scale_color_manual(values = c(GreenBay, Muskegon, GrandHaven, MichiganCity, Milwaukee, Naubinway, SuttonsBay))+
  theme_classic()+
  theme(text=element_text(color = "black", size = 15),axis.text=element_text(color = "black", size = 15))+
  theme(legend.position=c(0.7,0.85))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme(axis.ticks.length=unit(0.25, "cm"))+
  theme(axis.ticks=element_line(size=0.6,colour="black"))+
  theme(axis.line=element_line(size=0.6,colour="black"))
