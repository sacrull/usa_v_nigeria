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

#rarefied samples
rare_its <- rarefy_even_depth(otu_table(its.dat), rngseed = TRUE, replace = FALSE)
data_rare_its = data.frame(otu_table(rare_its)) # create a separated file
ps.rar.its <- phyloseq(rare_its, tax_its, map) # create a phyloseq object

#candida albicans
ps.rar.candida.its <- subset_taxa(ps.rar.its, V8 == "Candida_albicans")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.c_albicans.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.c_albicans.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()

#candida tropicalis
ps.rar.candida.its <- subset_taxa(ps.rar.its, V8 == "Candida_tropicalis")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.c_tropicalis.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.c_tropicalis.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()

#Debaryomyces_prosopidis
ps.rar.candida.its <- subset_taxa(ps.rar.its, V8 == "Debaryomyces_prosopidis")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.health.d_prosopidis.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Tooth_Classification))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()
pdf("./asv_dis.geo.d_prosopidis.pdf")
ggplot(data,aes(x=factor(names),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = order_50, guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
dev.off()

