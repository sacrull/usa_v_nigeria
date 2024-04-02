#install packages
#install.packages("ggplot2")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Rhdf5lib", force =TRUE )
# BiocManager::install("phyloseq")
# BiocManager::install("microbiome")
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
# remotes::install_github("gauravsk/ranacapa")
# remotes::install_github("bryandmartin/CORNCOB", force = TRUE)
# install.packages("reshape2")
# install.packages("tidyverse")
# install.packages("ggpubr")
# remotes::install_github("phytomosaic/ecole")

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
library(smplot2)

#set seed
set.seed(12349)
setwd("~/usa_nigeria/diversity")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

#its
#taxonomy barchart in america vs nigeria
rel.abund <- microbiome::transform(its.dat, 'compositional') # get relative abundance
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[2]) # collapse 
data <- psmelt(glom) # create dataframe from phyloseq object
data$V3 <- as.character(data$V3) # convert to character
data$V3[data$Abundance < 0.05] <- "< 5% abund"
avs <- plyr::ddply(data, ~V3, function(x) c(mean=mean(x$Abundance)))
avs
grouped <- data %>% group_by(Geog_loc, V3) %>% summarize(Abundance = mean(Abundance))
# get correct order of group factors
pdf("./grouped_tax_barchart.its.pdf")
ggplot(grouped, aes(fill=V3, y=Abundance, x=Geog_loc)) + geom_bar(position="fill", stat="identity") + theme_minimal()
dev.off()

#normalization
#clr transform
its.dat.clr <- microbiome::transform(its.dat, transform="clr", target="OTU")
#betadiversity
#color pallette
geo_color <- c("#8213A0", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

#cap plot by location
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.its.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean") # check signficane
#adonis beta diveristy
bd_clr_its <- phyloseq::distance(its.dat.clr, method = "euclidean")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(its.dat.clr))
adonis2(bd_clr_its ~ Geog_loc, data = sampledf, permutations=999)

#cap plot by tooth health
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~Tooth_Classification)
pdf("./bdiv_cap.tooth_health.its.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Tooth_Classification") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Tooth_Classification, method="euclidean") # check signficane
#adonis2 check sig
adonis2(bd_clr_its ~ Tooth_Classification, data = sampledf, permutations=999)

#seperate health into two plots: American and Nigeria
#usa
hlth_loco <- subset_samples(its.dat.clr, Geog_loc == "USA")
ordcap <- ordinate(hlth_loco, "CAP", "euclidean", ~Tooth_Classification)
pdf("./bdiv_cap.tooth_health.USA.its.pdf")
plot_ordination(hlth_loco, ordcap, "samples", color="Tooth_Classification") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()
permanova_pairwise(otu_table(hlth_loco), grp=sample_data(hlth_loco)$Tooth_Classification, method="euclidean") # check signficane
#nigeria
hlth_loco <- subset_samples(its.dat.clr, Geog_loc == "Nigeria" & Sample_Name != "DM00055V1PQ75") #sample was removed for being outlier
ordcap <- ordinate(hlth_loco, "CAP", "euclidean", ~Tooth_Classification)
pdf("./bdiv_cap.tooth_health.nigeria.its.pdf")
plot_ordination(hlth_loco, ordcap, "samples", color="Tooth_Classification") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()
permanova_pairwise(otu_table(hlth_loco), grp=sample_data(hlth_loco)$Tooth_Classification, method="euclidean") # check signficane

#beta dispersion
# first pull sample data from phyloseq object
metadata <- as(sample_data(its.dat.clr), "data.frame")
# calculate aitchison distance (from CLR transformed data)
clr.dist <- dist(otu_table(its.dat.clr), method="euclidean")
#geo locations
dispr <- vegan::betadisper(clr.dist, phyloseq::sample_data(its.dat.clr)$Geog_loc)
print("Beta disperson HIV status")
dispr
permutest(dispr)
pdf("./bdisp.geo_loc.its.pdf")
	boxplot(dispr)
dev.off()
#tooth health
dispr <- vegan::betadisper(clr.dist, phyloseq::sample_data(its.dat.clr)$Tooth_Classification)
print("Beta disperson HIV status")
dispr
permutest(dispr)
pdf("./bdisp.tooth_health.its.pdf")
	boxplot(dispr)
dev.off()
#permanova
metadata <- as(sample_data(its.dat.clr), "data.frame")
adonis2(clr.dist ~ Tooth_Classification * Geog_loc, data=metadata)


#alpha diversity
#nontransformed samples
#geo location
pdf("./adiv.geo_loc.its.pdf")
plot_richness(its.dat, measures=c("Observed", "Shannon"), x="Geog_loc") + 
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
#tooth health
pdf("./adiv.tooth_health.its.pdf")
plot_richness(its.dat, measures=c("Observed", "Shannon"), x="Tooth_Classification") + 
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

#rarefied samples
rare_its <- rarefy_even_depth(otu_table(its.dat), rngseed = TRUE, replace = FALSE)
otu.rare.its = data.frame(otu_table(rare_its)) # create a separated file
its.rare <- phyloseq(rare_its, tax_its, refseq_its, map) # create a phyloseq object

pdf("./adivrare.geo_loc.its.pdf")
plot_richness(its.rare, measures=c("Observed", "Shannon"), x="Geog_loc") + 
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
#tooth health
pdf("./adivrare.tooth_health.its.pdf")
plot_richness(its.rare, measures=c("Observed", "Shannon"), x="Tooth_Classification") + 
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
#rarefied
rare_its <- rarefy_even_depth(otu_table(its.dat), rngseed = TRUE, replace = FALSE)
data_rare_its = data.frame(otu_table(rare_its)) # create a separated file
ps.rar.its <- phyloseq(rare_its, tax_its, map) # create a phyloseq object
#geo location
pdf("./adiv_rare.geo_loc.its.pdf")
plot_richness(ps.rar.its, measures=c("Observed", "Shannon"), x="Geog_loc") + 
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
#tooth health
pdf("./adiv_rare.tooth_health.its.pdf")
plot_richness(ps.rar.rpoc, measures=c("Observed", "Shannon"), x="Tooth_Classification") + 
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
#rarefaction
options(warn=-1) # suppress warnings
p <- ggrare(its.dat, step = 1000, color = "Geog_loc", se = TRUE)
p <- p + facet_wrap(~Tooth_Classification)
pdf("./rarefaction_plots.its.pdf")
p + theme_minimal()
dev.off()
options(warn=0) # back on


# proportion of ASV reads found in geo location
top50 <-paste0("ASV",1:54)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object

pdf("./asv_dis.its.pdf")
ggplot(data,aes(x=factor(OTU),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = top50,guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

#upset plot
merged <- merge(otu_table(its.dat), sample_data(its.dat), by="row.names")
n <- ncol(otu_table(its.dat)) + 1
# by sex
agg <- aggregate(merged[,2:n], by=list(merged$Geog_loc), FUN=sum) 
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table 
agg[agg>1] <- 1
agg <- data.frame(t(agg[,-1]))
pdf("./upset.geo.its.pdf")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs")
dev.off()


#rpoc
#taxonomy barchart in america vs nigeria
rel.abund <- microbiome::transform(rpoc.dat, 'compositional') # get relative abundance
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[2]) # collapse 
data <- psmelt(glom) # create dataframe from phyloseq object
data$V3 <- as.character(data$V3) # convert to character
data$V3[data$Abundance < 0.05] <- "< 5% abund"
avs <- plyr::ddply(data, ~V3, function(x) c(mean=mean(x$Abundance)))
avs
grouped <- data %>% group_by(Geog_loc, V3) %>% summarize(Abundance = mean(Abundance))
# get correct order of group factors
pdf("./grouped_rpoc_tax_barchart.pdf")
ggplot(grouped, aes(fill=V3, y=Abundance, x=Geog_loc)) + geom_bar(position="fill", stat="identity") + theme_minimal()
dev.off()

#normalization
#clr transform
rpoc.dat.clr <- microbiome::transform(rpoc.dat, transform="clr", target="OTU")
#betadiversity
#color pallette
geo_color <- c("#8213A0", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

#cap plot by location
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo_loc.rpoc.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean") # check signficane
#euc by location
ordeuc <- ordinate(rpoc.dat.clr, "PCoA", "euclidean")
pdf("./bdiv_euc.geo_loc.rpoc.pdf")
plot_ordination(rpoc.dat.clr, ordeuc, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
#cap plot by tooth health
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Tooth_Classification)
pdf("./bdiv_cap.tooth_health.rpoc.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Tooth_Classification", shape ="Geog_loc") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()

metadata <- as(sample_data(rpoc.dat.clr), "data.frame")
adonis2(clr.dist ~ Tooth_Classification * Geog_loc, data=metadata)
#seperate health into two plots: American and Nigeria
#usa
hlth_loco <- subset_samples(rpoc.dat.clr, Geog_loc == "USA")
ordcap <- ordinate(hlth_loco, "CAP", "euclidean", ~Tooth_Classification)
pdf("./bdiv_cap.tooth_health.USA.rpoc.pdf")
plot_ordination(hlth_loco, ordcap, "samples", color="Tooth_Classification") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()
permanova_pairwise(otu_table(hlth_loco), grp=sample_data(hlth_loco)$Tooth_Classification, method="euclidean") # check signficane
#nigeria
hlth_loco <- subset_samples(rpoc.dat.clr, Geog_loc == "Nigeria") #sample was removed for being outlier
ordcap <- ordinate(hlth_loco, "CAP", "euclidean", ~Tooth_Classification)
pdf("./bdiv_cap.tooth_health.nigeria.rpoc.pdf")
plot_ordination(hlth_loco, ordcap, "samples", color="Tooth_Classification") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()
permanova_pairwise(otu_table(hlth_loco), grp=sample_data(hlth_loco)$Tooth_Classification, method="euclidean") # check signficane

#beta dispersion
# first pull sample data from phyloseq object
metadata <- as(sample_data(rpoc.dat.clr), "data.frame")
# calculate aitchison distance (from CLR transformed data)
clr.dist <- dist(otu_table(rpoc.dat.clr), method="euclidean")
#geo locations
dispr <- vegan::betadisper(clr.dist, phyloseq::sample_data(rpoc.dat.clr)$Geog_loc)
print("Beta disperson HIV status")
dispr
permutest(dispr)
pdf("./bdisp.geo_loc.rpoc.pdf")
    boxplot(dispr)
dev.off()
#tooth health
dispr <- vegan::betadisper(clr.dist, phyloseq::sample_data(rpoc.dat.clr)$Tooth_Classification)
print("Beta disperson HIV status")
dispr
permutest(dispr)
pdf("./bdisp.tooth_health.rpoc.pdf")
    boxplot(dispr)
dev.off()
#permanova
metadata <- as(sample_data(rpoc.dat.clr), "data.frame")
adonis2(clr.dist ~ Tooth_Classification * Geog_loc, data=metadata)


#alpha diversity
#nontransformed samples
#geo location
pdf("./adiv.geo_loc.rpoc.pdf")
plot_richness(rpoc.dat, measures=c("Observed", "Shannon"), x="Geog_loc") + 
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
#tooth health
pdf("./adiv.tooth_health.rpoc.pdf")
plot_richness(rpoc.dat, measures=c("Observed", "Shannon"), x="Tooth_Classification") + 
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
#rarefied
rare_rpoc <- rarefy_even_depth(otu_table(rpoc.dat), rngseed = TRUE, replace = FALSE)
data_rare_rpoc = data.frame(otu_table(rare_rpoc)) # create a separated file
ps.rar.rpoc <- phyloseq(rare_rpoc, tax_rpoc, map) # create a phyloseq object
#geo location
pdf("./adiv_rare.geo_loc.rpoc.pdf")
plot_richness(ps.rar.rpoc, measures=c("Observed", "Shannon"), x="Geog_loc") + 
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
#tooth health
pdf("./adiv_rare.tooth_health.rpoc.pdf")
plot_richness(ps.rar.rpoc, measures=c("Observed", "Shannon"), x="Tooth_Classification") + 
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
#rarefaction
options(warn=-1) # suppress warnings
p <- ggrare(rpoc.dat, step = 1000, color = "Geog_loc", se = TRUE)
p <- p + facet_wrap(~Tooth_Classification)
pdf("./rarefaction_plots.rpoc.pdf")
p + theme_minimal()
dev.off()
options(warn=0) # back on

# proportion of ASV reads found in geo location
top50 <-paste0("rASV",1:51)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object

pdf("./asv_dis.rpoc.pdf")
ggplot(data,aes(x=factor(OTU),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = top50,guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

#upset plot
merged <- merge(otu_table(rpoc.dat), sample_data(rpoc.dat), by="row.names")
n <- ncol(otu_table(rpoc.dat)) + 1
# by sex
agg <- aggregate(merged[,2:n], by=list(merged$Geog_loc), FUN=sum) 
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table 
agg[agg>1] <- 1
agg <- data.frame(t(agg[,-1]))
pdf("./upset.geo.rpoc.pdf")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs")
dev.off()


#correlation alpha diversity of bact and fungi
alpha.rpoc <- as.data.frame(estimate_richness(rpoc.dat, split = TRUE, measures = "Shannon")) #caculate shannon diversity
alpha.its <- as.data.frame(estimate_richness(its.dat, split = TRUE, measures = "Shannon"))
alpha.rpoc$names <- rownames(alpha.rpoc)
alpha.its$names <- rownames(alpha.its)
colnames(alpha.rpoc)[1] ="rpoC"
colnames(alpha.its)[1] ="its"
alpha_div <- inner_join(alpha.rpoc, alpha.its, by=c('names'='names'))
meta <- select(metadata, Sample_Name, Tooth_Classification, Geog_loc)
colnames(meta)[1] ="names"
alpha_div <- inner_join(alpha_div, meta, by=c('names'='names'))

pdf("test.pdf")
 ggscatter(alpha_div, x = "rpoC", y = "its",
   color = "Geog_loc",
   add = "reg.line", conf.int = TRUE)+
 stat_cor(aes(color = Geog_loc), label.x = 3)
dev.off()



#messing around
ps_50 <- subset_taxa(rpoc.pd, V8=="Leptotrichia_sp._oral_taxon_212")
data <- psmelt(ps_50) # create dataframe from phyloseq object

pdf("./asv_dis.test.pdf")
ggplot(data,aes(x=factor(OTU),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()