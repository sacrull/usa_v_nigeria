#load libraries
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

#set seed
set.seed(12349)
setwd("~/usa_nigeria/hcf_comp")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

#clr transform
its.hcf <- subset_samples(its.dat, Tooth_Classification == "H-CF")
its.dat.clr <- microbiome::transform(its.hcf, transform="clr", target="OTU")

#cap plot by location
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.its_hcf.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean") # check signficane

#alpha diversity
#geo location
pdf("./adiv.geo_loc.its_hcf.pdf")
plot_richness(its.hcf, measures=c("Observed", "Shannon"), x="Geog_loc") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()


#clr transform
rpoc.hcf <- subset_samples(rpoc.dat, Tooth_Classification == "H-CF")
rpoc.dat.clr <- microbiome::transform(rpoc.hcf, transform="clr", target="OTU")

#cap plot by location
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.rpoc_hcf.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean") # check signficane

#alpha diversity
#geo location
pdf("./adiv.geo_loc.rpoc_hcf.pdf")
plot_richness(rpoc.hcf, measures=c("Observed", "Shannon"), x="Geog_loc") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()

#DCD
its.dcd <- subset_samples(its.dat, Tooth_Classification == "D-CD")
its.dat.clr <- microbiome::transform(its.dcd, transform="clr", target="OTU")

#cap plot by location
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.its_dcd.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean") # check signficane

#alpha diversity
#geo location
pdf("./adiv.geo_loc.its_dcd.pdf")
plot_richness(its.dcd, measures=c("Observed", "Shannon"), x="Geog_loc") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()


#clr transform
rpoc.dcd <- subset_samples(rpoc.dat, Tooth_Classification == "D-CD")
rpoc.dat.clr <- microbiome::transform(rpoc.dcd, transform="clr", target="OTU")

#cap plot by location
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.rpoc_dcd.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean") # check signficane

#alpha diversity
#geo location
pdf("./adiv.geo_loc.rpoc_dcd.pdf")
plot_richness(rpoc.dcd, measures=c("Observed", "Shannon"), x="Geog_loc") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()

library(grid)
library(coda4microbiome)
library(tidyverse)
#bal for HCF
# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(its.hcf, "V8")
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
dif2 <- dim(dat)[2] - 4
y <- factor(dat[,dif2]) #geog lcaotion
z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound
geo_its <- coda_glmnet(x=x,y=y) #running it
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
pdf("./geo.its_hcf.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(its.dcd, "V8")
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
dif2 <- dim(dat)[2] - 4
y <- factor(dat[,dif2]) #geog lcaotion
z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound
geo_its <- coda_glmnet(x=x,y=y) #running it
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
pdf("./geo.its_dcd.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()

#rpoc
# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(rpoc.hcf, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
#get taxonomy of ASV
taxa = as(tax_table(rpoc.dat), "matrix")
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
y <- factor(dat[,dif2]) #geog lcaotion
z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound
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
pdf("./geo.rpoc_hcf.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(rpoc.dcd, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
#glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)
#get taxonomy of ASV
taxa = as(tax_table(rpoc.dat), "matrix")
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
y <- factor(dat[,dif2]) #geog lcaotion
z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound
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
pdf("./geo.rpoc_dcd.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()