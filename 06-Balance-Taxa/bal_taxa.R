#devtools::install_github(repo = "malucalle/selbal")

library(selbal)
library(grid)
library(coda4microbiome)
library(tidyverse)
library(phyloseq)


setwd("~/usa_nigeria/bal_tax")
load("~/usa_nigeria/phyloseq_obj/ps.RData")
set.seed(12349)

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(its.dat, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
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

pdf("./geo.its.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()


# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(rpoc.dat, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
# pull data
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

pdf("./geo.rpoc.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()

