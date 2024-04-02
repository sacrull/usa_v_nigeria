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
setwd("~/usa_nigeria/all_cats/diversity")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

#its
its.dat <- prune_samples(sample_sums(its.pd) > 2000, its.pd) #remove less than 2000 reads
its.dat <- subset_taxa(its.dat, V2!="Eukaryote_unknown" & V3!="Fungi_unknown") #remove unwanted taxa

#rpoc
rpoc.dat <- prune_samples(sample_sums(rpoc.pd) > 4000, rpoc.pd)

#its
its.dat.clr <- microbiome::transform(its.dat, transform="clr", target="OTU")
#betadiversity
#color pallette
geo_color <- c("#8213A0", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

#cap plot by location
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~hiv_status)
pdf("./bdiv_cap.hiv_status.its.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="hiv_status") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
sig <- subset_samples(its.dat.clr, tooth_health == "H")
permanova_pairwise(otu_table(sig), grp=sample_data(sig)$hiv_status, method="euclidean") # check signficane
sig <- subset_samples(its.dat.clr, tooth_health == "D")
permanova_pairwise(otu_table(sig), grp=sample_data(sig)$hiv_status, method="euclidean") # check signficane
#unsuprivised beta diveristy
ordeuc <- ordinate(its.dat.clr, "PCoA", "euclidean", ~hiv_status)
pdf("./bdiv_euc.hiv.its_mouth.pdf")
plot_ordination(its.dat.clr, ordeuc, "samples", color="hiv_status") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~mouth_health)
dev.off()

#upset plot
merged <- merge(otu_table(its.dat), sample_data(its.dat), by="row.names")
n <- ncol(otu_table(its.dat)) + 1
# by sex
agg <- aggregate(merged[,2:n], by=list(merged$hiv_status), FUN=sum) 
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table 
agg[agg>1] <- 1
agg <- data.frame(t(agg[,-1]))
pdf("./upset.hiv.its.pdf")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs")
dev.off()

#alpha diversity
pdf("./adiv.hiv_stat.its.pdf")
plot_richness(its.dat, measures=c("Observed", "Shannon"), x="hiv_status") + 
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

#rpoc
rpoc.dat.clr <- microbiome::transform(rpoc.dat, transform="clr", target="OTU")
#betadiversity
#color pallette
geo_color <- c("#8213A0", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

#cap plot by location
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~hiv_status)
pdf("./bdiv_cap.hiv_status.rpoc.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="hiv_status") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$hiv_status, method="euclidean") # check signficane
#unsuprivised beta diveristy
ordeuc <- ordinate(rpoc.dat.clr, "PCoA", "euclidean", ~hiv_status)
pdf("./bdiv_euc.hiv.rpoc_mouth.pdf")
plot_ordination(rpoc.dat.clr, ordeuc, "samples", color="hiv_status") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)+
    facet_wrap(~mouth_health)
dev.off()

#upset plot
merged <- merge(otu_table(rpoc.dat), sample_data(rpoc.dat), by="row.names")
n <- ncol(otu_table(rpoc.dat)) + 1
# by sex
agg <- aggregate(merged[,2:n], by=list(merged$hiv_status), FUN=sum) 
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table 
agg[agg>1] <- 1
agg <- data.frame(t(agg[,-1]))
pdf("./upset.hiv.rpoc.pdf")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs")
dev.off()

#alpha diversity
pdf("./adiv.hiv_stat.rpoc.pdf")
plot_richness(rpoc.dat, measures=c("Observed", "Shannon"), x="hiv_status") + 
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