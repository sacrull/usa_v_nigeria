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
setwd("~/usa_nigeria/bal_tax/all_samples")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

its.dat <- prune_samples(sample_sums(its.pd) > 2000, its.pd) #remove less than 2000 reads
its.dat <- subset_samples(its.dat , hiv_status == "HUU" | hiv_status == "NAN")
its.dat <- subset_taxa(its.dat, V2!="Eukaryote_unknown" & V3!="Fungi_unknown")
#usa
# collapse data to roughly species level to minimize high sparsity
its.geo <- subset_samples(its.dat, Geog_loc == "USA")
its.cd <- subset_samples(its.geo, tooth_health == "H" | tooth_health == "D")
glom <- tax_glom(its.cd, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
#get taxonomy of ASV
taxa = as(tax_table(its.pd), "matrix")
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
dif2 <- dim(dat)[2] - 4
y <- factor(dat[,dif2]) #H or D
#running coda
geo_rpoc <- coda_glmnet(x=x,y=y) 
sum(geo_rpoc$`log-contrast coefficients`)
#positive taxa
coef<-geo_rpoc$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_rpoc$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_rpoc$taxa.name[negatives[on]]
pdf("./its.H_D.usa.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()
#nigeria
# collapse data to roughly species level to minimize high sparsity
its.geo <- subset_samples(its.dat, Geog_loc == "Nigeria" & hiv_status == "HUU")
its.cd <- subset_samples(its.geo, tooth_health == "H" | tooth_health == "D")
glom <- tax_glom(its.cd, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
#get taxonomy of ASV
taxa = as(tax_table(its.pd), "matrix")
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
dif2 <- dim(dat)[2] - 4
y <- factor(dat[,dif2]) #H or D
#running coda
geo_rpoc <- coda_glmnet(x=x,y=y) 
sum(geo_rpoc$`log-contrast coefficients`)
#positive taxa
coef<-geo_rpoc$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_rpoc$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_rpoc$taxa.name[negatives[on]]
pdf("./its.H_D.nigeria.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()
#H USA vs Nigeria
# collapse data to roughly species level to minimize high sparsity
its.geo <- subset_samples(its.dat, hiv_status == "NAN" | hiv_status == "HUU")
its.cd <- subset_samples(its.geo, tooth_health == "H")
glom <- tax_glom(its.cd, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
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
geo_rpoc <- coda_glmnet(x=x,y=y) 
sum(geo_rpoc$`log-contrast coefficients`)
#positive taxa
coef<-geo_rpoc$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_rpoc$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_rpoc$taxa.name[negatives[on]]
pdf("./its.usa_nigeria.H.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()
#USA vs Nigeria D
# collapse data to roughly species level to minimize high sparsity
its.geo <- subset_samples(its.dat, hiv_status == "NAN" | hiv_status == "HUU")
its.cd <- subset_samples(its.geo, tooth_health == "D")
glom <- tax_glom(its.cd, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
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
geo_rpoc <- coda_glmnet(x=x,y=y) 
sum(geo_rpoc$`log-contrast coefficients`)
#positive taxa
coef<-geo_rpoc$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_rpoc$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_rpoc$taxa.name[negatives[on]]
pdf("./its.usa_nigeria.D.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()