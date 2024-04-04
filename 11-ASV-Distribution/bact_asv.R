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
setwd("~/usa_nigeria/candida")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

rare_rpoc <- rarefy_even_depth(otu_table(rpoc.dat), rngseed = TRUE, replace = FALSE)
data_rare_rpoc = data.frame(otu_table(rare_rpoc)) # create a separated file
ps.rar.rpoc <- phyloseq(rare_rpoc, tax_rpoc, map) # create a phyloseq object
#strep mutans
ps.rar.candida.rpoc <- subset_taxa(ps.rar.rpoc, V8 == "Streptococcus_mutans")
#more than 10 reads in 5% or more samples
#ps.rar.candida.rpoc <- filter_taxa(ps.rar.candida.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.rpoc, n=50)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.s_mutans.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.s_mutans.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()

#lepto 212
ps.rar.candida.rpoc <- subset_taxa(ps.rar.rpoc, V8 == "Leptotrichia_sp._oral_taxon_212")
#more than 10 reads in 5% or more samples
#ps.rar.candida.rpoc <- filter_taxa(ps.rar.candida.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.rpoc, n=50)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.lepto_212.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.lepto_212.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()

#lepto 215
ps.rar.candida.rpoc <- subset_taxa(ps.rar.rpoc, V8 == "Leptotrichia_sp._oral_taxon_215")
#more than 10 reads in 5% or more samples
#ps.rar.candida.rpoc <- filter_taxa(ps.rar.candida.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.rpoc, n=50)
ps_50 <- subset_taxa(ps.rar.rpoc, rownames(tax_table(ps.rar.rpoc)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.lepto_215.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.lepto_215.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()


#trepnema denticola
ps.rar.candida.rpoc <- subset_taxa(rpoc.dat, V8 == "Treponema_denticola")
#more than 10 reads in 5% or more samples
#ps.rar.candida.rpoc <- filter_taxa(ps.rar.candida.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.rpoc, n=50)
ps_50 <- subset_taxa(rpoc.dat, rownames(tax_table(rpoc.dat)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.treponema_denticola.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.treponema_denticola.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()

#Porphyromonas_gingivalis
ps.rar.candida.rpoc <- subset_taxa(rpoc.dat, V8 == "Porphyromonas_gingivalis")
#more than 10 reads in 5% or more samples
#ps.rar.candida.rpoc <- filter_taxa(ps.rar.candida.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.rpoc, n=50)
ps_50 <- subset_taxa(rpoc.dat, rownames(tax_table(rpoc.dat)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.porphyromonas_gingivalis.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.porphyromonas_gingivalis.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()

#Tannerella forsythia
ps.rar.candida.rpoc <- subset_taxa(rpoc.dat, V8 == "Tannerella_forsythia")
#more than 10 reads in 5% or more samples
#ps.rar.candida.rpoc <- filter_taxa(ps.rar.candida.rpoc, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.rpoc, n=50)
ps_50 <- subset_taxa(rpoc.dat, rownames(tax_table(rpoc.dat)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.tannerella_forsythia.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.tannerella_forsythia.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()