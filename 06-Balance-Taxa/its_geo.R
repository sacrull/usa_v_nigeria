#its
library(ggplot2, verbose=F)
library(phyloseq, verbose=F)
library(ape, verbose=F)
library(metagMisc, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(vegan, verbose=F)
library(ranacapa, verbose=F)
library(microbiome, verbose=F)
library(corncob, verbose=F)
library(magrittr, verbose=F)
library(ggpubr, verbose=F)
library(ecole, verbose=F)
library(UpSetR)
library(stringr)
library(reshape2)
library(microbiome)
library(grid)
library(coda4microbiome)
library(tidyverse)

#set seed
set.seed(12349)
setwd("~/usa_nigeria/bal_tax")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

#tooth health
# collapse data to roughly species level to minimize high sparsity
its.geo <- subset_samples(its.dat, tooth_health == "H")
glom <- tax_glom(its.geo, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
#get taxonomy of ASV
taxa = as(tax_table(its.dat), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
#rename to have V8 level name
dat <- as.data.frame(dat)
dat <- dat %>% rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))
rownames(dat) <- dat$V8
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names
# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,1:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# response variable
dif2 <- dim(dat)[2] - 5
y <- factor(dat[,dif2]) #usa or nigeria
#running coda
geo_its <- coda_glmnet(x=x,y=y) 
sum(geo_its$`log-contrast coefficients`)
#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]
pdf("./its.h.geo.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()

#d
# collapse data to roughly species level to minimize high sparsity
its.geo <- subset_samples(its.dat, tooth_health == "D")
glom <- tax_glom(its.geo, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
#get taxonomy of ASV
taxa = as(tax_table(its.dat), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
#rename to have V8 level name
dat <- as.data.frame(dat)
dat <- dat %>% rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))
rownames(dat) <- dat$V8
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names
# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,1:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# response variable
dif2 <- dim(dat)[2] - 5
y <- factor(dat[,dif2]) #usa or nigeria
#running coda
geo_its <- coda_glmnet(x=x,y=y) 
sum(geo_its$`log-contrast coefficients`)
#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]
pdf("./its.d.geo.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()