#rpoc
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
rpoc.hcf <- subset_samples(rpoc.dat, Tooth_Classification == "H-CF")
rpoc.dcd <- subset_samples(rpoc.dat, Tooth_Classification == "D-CD")
rpoc.dat.clr <- microbiome::transform(rpoc.dat, transform="clr", target="OTU")

#cap plot by location
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Geog_loc+Tooth_Classification)
pdf("./bdiv_cap.geo_loc.rpoc_mouth.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~mouth_health)
dev.off()
pdf("./bdiv_cap.geo_loc.rpoc_tooth.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~Tooth_Classification)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean") # check signficane

ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.rpoc_tooth_health.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~tooth_health)
dev.off()
rpoc.usa <- subset_samples(rpoc.dat.clr, Geog_loc == "USA")
permanova_pairwise(otu_table(rpoc.usa), grp=sample_data(rpoc.usa)$tooth_health, method="euclidean") # check signficane
rpoc.nigeria <- subset_samples(rpoc.dat.clr, Geog_loc == "Nigeria")
permanova_pairwise(otu_table(rpoc.nigeria), grp=sample_data(rpoc.nigeria)$tooth_health, method="euclidean") # check signficane
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~mouth_health)
pdf("./bdiv_cap.geo_loc.rpoc_mouth_health.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="mouth_health") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~Geog_loc)
dev.off()

#unsupervised beta diversity
ordeuc <- ordinate(rpoc.dat.clr, "PCoA", "euclidean")
pdf("./bdiv_euc.geo_loc.rpoc_mouth.pdf")
plot_ordination(rpoc.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~mouth_health)
dev.off()
pdf("./bdiv_euc.geo_loc.rpoc_tooth_health.pdf")
plot_ordination(rpoc.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~tooth_health)
dev.off()
pdf("./bdiv_euc.geo_loc.rpoc_tooth.pdf")
plot_ordination(rpoc.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~Tooth_Classification)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean") # check signficane
ordeuc <- ordinate(its.dat.clr, "PCoA", "euclidean")
pdf("./bdiv_euc.geo_loc.rpoc_tooth_health.pdf")
plot_ordination(rpoc.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~tooth_health)
dev.off()
rpoc.tooth <- subset_samples(rpoc.dat.clr, tooth_health == "H")
permanova_pairwise(otu_table(rpoc.tooth), grp=sample_data(rpoc.tooth)$Geog_loc, method="euclidean") # check signficane
rpoc.tooth <- subset_samples(rpoc.dat.clr, tooth_health == "D")
permanova_pairwise(otu_table(rpoc.tooth), grp=sample_data(rpoc.tooth)$Geog_loc, method="euclidean") # check signficane



#comparing health categories
rpoc.hcf <- subset_samples(rpoc.dat, Tooth_Classification == "H-CF")
rpoc.dat.clr <- microbiome::transform(rpoc.hcf, transform="clr", target="OTU")
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean", padj = 'BH') # check signficane
rpoc.dcd <- subset_samples(rpoc.dat, Tooth_Classification == "D-CD")
rpoc.dat.clr <- microbiome::transform(rpoc.dcd, transform="clr", target="OTU")
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean", padj = 'BH') # check signficane

#alpha diversity
#geo location
pdf("./adiv.geo_loc.rpoc_tooth.pdf")
plot_richness(rpoc.dat, measures=c("Observed", "Shannon"), x="Geog_loc") + 
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

pdf("./adiv.geo_loc.rpoc_tooth.pdf")
plot_richness(rpoc.dat, measures=c("Observed", "Shannon"), x="Geog_loc") + 
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


pdf("./adiv.rpoc.health.pdf")
plot_richness(rpoc.dat, measures=c("Observed", "Shannon"), x="tooth_health") + 
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

pdf("./adiv.rpoc.oral.pdf")
plot_richness(rpoc.dat, measures=c("Observed", "Shannon"), x="Tooth_Classification") + 
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












rare_rpoc <- rarefy_even_depth(otu_table(rpoc.dat), rngseed = TRUE, replace = FALSE)
data_rare_rpoc = data.frame(otu_table(rare_rpoc)) # create a separated file
ps.rar.rpoc <- phyloseq(rare_rpoc, tax_rpoc, map) # create a phyloseq object

#usa
ps.rar.geo.rpoc <- subset_samples(ps.rar.rpoc, Geog_loc == "USA") # create a phyloseq object
ps.rar.geo.rpoc <- filter_taxa(ps.rar.geo.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.rpoc, n=50)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names

pdf("./asv_dis.rpoc_mouth.usa.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
pdf("./asv_dis.rpoc_tooth.usa.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
#nigeria
ps.rar.geo.rpoc <- subset_samples(ps.rar.rpoc, Geog_loc == "Nigeria") # create a phyloseq object
ps.rar.geo.rpoc <- filter_taxa(ps.rar.geo.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.rpoc, n=50)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names

pdf("./asv_dis.rpoc_mouth.nigeria.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
pdf("./asv_dis.rpoc_tooth.nigeria.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()

#glommed
glom <- tax_glom(ps.rar.rpoc, "V8")
ps.rar.geo.rpoc <- subset_samples(glom, Geog_loc == "USA") # create a phyloseq object
ps.rar.geo.rpoc <- filter_taxa(ps.rar.geo.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.rpoc, n=50)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, V8, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$V8

pdf("./tax_dis.rpoc_mouth.usa.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
pdf("./tax_dis.rpoc_tooth.usa.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()

ps.rar.geo.rpoc <- subset_samples(glom, Geog_loc == "Nigeria") # create a phyloseq object
ps.rar.geo.rpoc <- filter_taxa(ps.rar.geo.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.geo.rpoc, n=50)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, V8, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$V8

pdf("./tax_dis.rpoc_mouth.nigeria.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(tooth_health))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits= order_50, guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
pdf("./tax_dis.rpoc_tooth.nigeria.pdf")
ggplot(data,aes(x=factor(V8),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()

















#health relationship 
perma <- permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Tooth_Classification, method="euclidean", padj = "fdr") # check signficane
#changing formated into matrix
perma$pairs <- gsub(' vs ', '_', perma$pairs) # perparing to split column
perma[c('cat1', 'cat2')] <- str_split_fixed(perma$pairs, '_', 2) #splitting vs column
long_df <- select(perma, cat1, cat2, R2) #getting cat1 and cat2 for health and R2 for relationship
long_df2 <- select(perma, cat2, cat1, R2)
names(long_df2)[1] <- "cat1"
names(long_df2)[2] <- "cat2"
same_cats <- data.frame(cat1=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'), 
                               cat2=c('H-CF', 'H-CE', 'H-CD', 'E-CE', 'E-CD', 'D-CD'),
                               R2=c(1, 1, 1, 1, 1, 1), 
                               stringsAsFactors=FALSE) #finishing up the columns
long_df1 <- rbind(long_df, same_cats, long_df2) #combining dataframes
long_df3 <- dcast(long_df1, cat1 ~ cat2) #geting matrix of R2 values
mat <- as.matrix(data.frame(long_df3, row.names = 1))
my_nj <- ape::nj(mat) #making tree
write.tree(my_nj, file = "health_cat.rpoc.tre", append = FALSE,
                digrpoc = 10, tree.names = FALSE)
pdf("tooth_health_rel.rpoc.pdf")
plot(my_nj, "unrooted")
dev.off()
#by geography
rpoc.geo <- subset_samples(rpoc.dat, Geog_loc == "USA")
rpoc.geo.clr <- microbiome::transform(rpoc.geo, transform="clr", target="OTU")
perma <- permanova_pairwise(otu_table(rpoc.geo.clr), grp=sample_data(rpoc.geo.clr)$Tooth_Classification, method="euclidean", padj = "fdr")
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
write.tree(my_nj, file = "health_cat.rpoc.usa.tre", append = FALSE,
                digits = 10, tree.names = FALSE)
rpoc.geo <- subset_samples(rpoc.dat, Geog_loc == "Nigeria")
rpoc.geo.clr <- microbiome::transform(rpoc.geo, transform="clr", target="OTU")
perma <- permanova_pairwise(otu_table(rpoc.geo.clr), grp=sample_data(rpoc.geo.clr)$Tooth_Classification, method="euclidean", padj = "fdr")
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
write.tree(my_nj, file = "health_cat.rpoc.nigeria.tre", append = FALSE,
                digits = 10, tree.names = FALSE)


#metadata correlations
env <- sample_data(rpoc.dat.clr)
env.mat <- env[, c(4,6,7,8,9)]
all.sum.rpoc <- bioenv((otu_table(rpoc.dat.clr)), env.mat, index = "euclidean")

#balance of taxas
library(grid)
library(coda4microbiome)
library(tidyverse)
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

#mouth health
# collapse data to roughly species level to minimize high sparsity
rpoc.cd <- subset_samples(rpoc.dat, mouth_health == "CD" | mouth_health == "CF")
glom <- tax_glom(rpoc.cd, "V8")
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
dif2 <- dim(dat)[2] - 3
y <- factor(dat[,dif2]) #CF or CD
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
pdf("./rpoc_cd.pdf")
geo_rpoc$`signature plot`
geo_rpoc$`predictions plot`
dev.off()