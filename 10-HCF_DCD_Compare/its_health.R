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
library(stringr)
library(reshape2)
library(microbiome)


#set seed
set.seed(12349)
setwd("~/usa_nigeria/hcf_comp")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

#clr transform
its.hcf <- subset_samples(its.dat, Tooth_Classification == "H-CF")
its.dat.clr <- microbiome::transform(its.dat, transform="clr", target="OTU")

#cap plot by location
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.its_tooth.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~tooth_health)
dev.off()
pdf("./bdiv_cap.geo_loc.its_tooth.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~Tooth_Classification)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean") # check signficane
#unsupervised
ordeuc <- ordinate(its.dat.clr, "PCoA", "euclidean", ~Geog_loc+Tooth_Classification)
pdf("./bdiv_euc.geo_loc.its_mouth.pdf")
plot_ordination(its.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~mouth_health)
dev.off()
pdf("./bdiv_euc.geo_loc.its_tooth.pdf")
plot_ordination(its.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~Tooth_Classification)
dev.off()
ordeuc <- ordinate(its.dat.clr, "PCoA", "euclidean")
pdf("./bdiv_euc.geo_loc.its_tooth_health.pdf")
plot_ordination(its.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~tooth_health)
dev.off()
its.tooth <- subset_samples(its.dat.clr, tooth_health == "H")
permanova_pairwise(otu_table(its.tooth), grp=sample_data(its.tooth)$Geog_loc, method="euclidean") # check signficane
its.tooth <- subset_samples(its.dat.clr, tooth_health == "D")
permanova_pairwise(otu_table(its.tooth), grp=sample_data(its.tooth)$Geog_loc, method="euclidean") # check signficane

ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~tooth_health)
pdf("./bdiv_cap.geo_loc.its_tooth_health.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="tooth_health") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~Geog_loc)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$tooth_health, method="euclidean") # check signficane
its.usa <- subset_samples(its.dat.clr, Geog_loc == "USA")
permanova_pairwise(otu_table(its.usa), grp=sample_data(its.usa)$tooth_health, method="euclidean") # check signficane
permanova_pairwise(otu_table(its.usa), grp=sample_data(its.usa)$mouth_health, method="euclidean") # check signficane
its.nigeria <- subset_samples(its.dat.clr, Geog_loc == "Nigeria")
permanova_pairwise(otu_table(its.nigeria), grp=sample_data(its.nigeria)$tooth_health, method="euclidean") # check signficane
permanova_pairwise(otu_table(its.nigeria), grp=sample_data(its.nigeria)$mouth_health, method="euclidean") # check signficane

#alpha diversity
#geo location
pdf("./adiv.geo_loc.its_mouth.pdf")
plot_richness(its.dat, measures=c("Observed", "Shannon"), x="Geog_loc") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    facet_wrap(~tooth_health)+
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()
pdf("./adiv.geo_loc.its_tooth.pdf")
plot_richness(its.dat, measures=c("Observed", "Shannon"), x="Geog_loc") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    facet_wrap(~Tooth_Classification)+
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()


#comparing health categories
its.hcf <- subset_samples(its.dat, Tooth_Classification == "H-CF")
its.dat.clr <- microbiome::transform(its.hcf, transform="clr", target="OTU")
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean", padj = 'BH') # check signficane
its.dcd <- subset_samples(its.dat, Tooth_Classification == "D-CD")
its.dat.clr <- microbiome::transform(its.dcd, transform="clr", target="OTU")
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean", padj = 'BH') # check signficane


#alpha diversity
#geo location
pdf("./adiv.its.health.pdf")
plot_richness(its.dat, measures=c("Observed", "Shannon"), x="tooth_health") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25)+
    facet_wrap(~Geog_loc) +
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

pdf("./adiv.its.oral.pdf")
plot_richness(its.dat, measures=c("Observed", "Shannon"), x="Tooth_Classification") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25)+
    facet_wrap(~Geog_loc) +
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


#asv distrubtion in health
#rarefied samples
rare_its <- rarefy_even_depth(otu_table(its.dat), rngseed = TRUE, replace = FALSE)
data_rare_its = data.frame(otu_table(rare_its)) # create a separated file
ps.rar.its <- phyloseq(rare_its, tax_its, map) # create a phyloseq object

#usa
ps.rar.geo.its <- subset_samples(ps.rar.its, Geog_loc == "USA")
#more than 10 reads in 5% or more samples
ps.rar.geo.its <- filter_taxa(ps.rar.geo.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names

pdf("./asv_dis.its_mouth.usa.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
pdf("./asv_dis.its_tooth.usa.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
#nigeria
ps.rar.geo.its <- subset_samples(ps.rar.its, Geog_loc == "Nigeria") # create a phyloseq object
ps.rar.geo.its <- filter_taxa(ps.rar.geo.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names

pdf("./asv_dis.its_mouth.nigeria.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

pdf("./asv_dis.its_tooth.nigeria.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()


#glommed
glom <- tax_glom(ps.rar.its, "V8")
ps.rar.geo.its <- subset_samples(glom, Geog_loc == "USA") # create a phyloseq object
ps.rar.geo.its <- filter_taxa(ps.rar.geo.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, V8, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$V8

pdf("./tax_dis.its_mouth.usa.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
pdf("./tax_dis.its_tooth.usa.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()

ps.rar.geo.its <- subset_samples(glom, Geog_loc == "Nigeria") # create a phyloseq object
ps.rar.geo.its <- filter_taxa(ps.rar.geo.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, V8, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$V8

pdf("./tax_dis.its_mouth.nigeria.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
pdf("./tax_dis.its_tooth.nigeria.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()



#health relationship 
perma <- permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Tooth_Classification, method="euclidean", padj = "fdr") # check signficane
#changing formated into matrix
perma$pairs <- gsub(' vs ', '_', perma$pairs)
perma[c('cat1', 'cat2')] <- str_split_fixed(perma$pairs, '_', 2)
long_df <- select(perma, cat1, cat2, R2)
long_df2 <- select(perma, cat2, cat1, R2)
names(long_df2)[1] <- "cat1"
names(long_df2)[2] <- "cat2"
same_cats <- data.frame(cat1=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'), 
                               cat2=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'),
                               R2=c(1, 1, 1, 1, 1, 1), 
                               stringsAsFactors=FALSE)
long_df1 <- rbind(long_df, same_cats, long_df2)
long_df3 <- dcast(long_df1, cat1 ~ cat2) #geting matrix of R2 values
mat <- as.matrix(data.frame(long_df3, row.names = 1))
my_nj <- ape::nj(mat)
write.tree(my_nj, file = "health_cat.its.tre", append = FALSE,
                digits = 10, tree.names = FALSE)
pdf("tooth_health_rel.its.pdf")
plot(my_nj, "unrooted")
dev.off()
#by geography
its.geo <- subset_samples(its.dat, Geog_loc == "USA")
its.geo.clr <- microbiome::transform(its.geo, transform="clr", target="OTU")
perma <- permanova_pairwise(otu_table(its.geo.clr), grp=sample_data(its.geo.clr)$Tooth_Classification, method="euclidean", padj = "fdr")
#changing formated into matrix
perma$pairs <- gsub(' vs ', '_', perma$pairs)
perma[c('cat1', 'cat2')] <- str_split_fixed(perma$pairs, '_', 2)
long_df <- select(perma, cat1, cat2, R2)
long_df2 <- select(perma, cat2, cat1, R2)
names(long_df2)[1] <- "cat1"
names(long_df2)[2] <- "cat2"
same_cats <- data.frame(cat1=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'), 
                               cat2=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'),
                               R2=c(1, 1, 1, 1, 1, 1), 
                               stringsAsFactors=FALSE)
long_df1 <- rbind(long_df, same_cats, long_df2)
long_df3 <- dcast(long_df1, cat1 ~ cat2) #geting matrix of R2 values
mat <- as.matrix(data.frame(long_df3, row.names = 1))
my_nj <- ape::nj(mat)
write.tree(my_nj, file = "health_cat.its.usa.tre", append = FALSE,
                digits = 10, tree.names = FALSE)
its.geo <- subset_samples(its.dat, Geog_loc == "Nigeria")
its.geo.clr <- microbiome::transform(its.geo, transform="clr", target="OTU")
perma <- permanova_pairwise(otu_table(its.geo.clr), grp=sample_data(its.geo.clr)$Tooth_Classification, method="euclidean", padj = "fdr")
#changing formated into matrix
perma$pairs <- gsub(' vs ', '_', perma$pairs)
perma[c('cat1', 'cat2')] <- str_split_fixed(perma$pairs, '_', 2)
long_df <- select(perma, cat1, cat2, R2)
long_df2 <- select(perma, cat2, cat1, R2)
names(long_df2)[1] <- "cat1"
names(long_df2)[2] <- "cat2"
same_cats <- data.frame(cat1=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'), 
                               cat2=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'),
                               R2=c(1, 1, 1, 1, 1, 1), 
                               stringsAsFactors=FALSE)
long_df1 <- rbind(long_df, same_cats, long_df2)
long_df3 <- dcast(long_df1, cat1 ~ cat2) #geting matrix of R2 values
mat <- as.matrix(data.frame(long_df3, row.names = 1))
my_nj <- ape::nj(mat)
write.tree(my_nj, file = "health_cat.its.nigeria.tre", append = FALSE,
                digits = 10, tree.names = FALSE)


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

# collapse data to roughly species level to minimize high sparsity
its.dat.geo <- subset_samples(its.dat, Geog_loc == "Nigeria")
its.dat.geo <- subset_samples(its.dat.geo, Tooth_Classification == "H-CF" | Tooth_Classification == "D-CD")
glom <- tax_glom(its.dat.usa, "V8")
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
dif2 <- dim(dat)[2] - 6
y <- factor(dat[,dif2]) #helth status
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
pdf("./geo.its_dcd.nigeria.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
