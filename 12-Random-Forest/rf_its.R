#install.packages('randomForest')
#install.packages("caTools")
library(phyloseq)
library(randomForest)
library(caret)
library(e1071)
library(reshape2)
library(tidyverse)
library(dplyr)
library(caTools)
#set seed and load data
set.seed(12349)
setwd("~/usa_nigeria/random_forest")
load("~/usa_nigeria/phyloseq_obj/ps.RData")


###its
# format data
asv_tab <- read.table("~/usa_nigeria/its/sequence_table_its.merged.txt", sep="\t", header=T, row.names=1, stringsAsFactors=F, comment.char="")
metadata <- read.table("~/usa_nigeria/Nigeria_USA_meta.txt", sep="\t", header=T, row.names=1, stringsAsFactors=T, comment.char="")
metadata$Geog_loc <- factor(metadata$Geog_loc)

# get relative abudance
asv_tab_norm <- sweep(asv_tab, 2, colSums(asv_tab), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
# remove nas
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]
set.seed(151)

# by geo location
asv_tab_var$var <- metadata[rownames(asv_tab_var), "Geog_loc"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group

# plot important ASVs
pdf("./rf.geo.importance.its.pdf")
varImpPlot(rf.study_group)
dev.off()

# by health location
metadata$Tooth_Classification <- factor(metadata$Tooth_Classification)
asv_tab_var$var <- metadata[rownames(asv_tab_var), "Tooth_Classification"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group

# plot important ASVs
pdf("./rf.tooth.importance.its.pdf")
varImpPlot(rf.study_group)
dev.off()

#within location
# USA
asv_tab2 <- asv_tab[!grepl("DM00", row.names(asv_tab)),] #remove nigeria samples
asv_tab_norm <- sweep(asv_tab2, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata2 <- metadata[!grepl("DM00", row.names(metadata)),] #remove nigeria samples
metadata2$Tooth_Classification <- factor(metadata2$Tooth_Classification)
asv_tab_var$var <- metadata2[rownames(asv_tab_var), "Tooth_Classification"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.tooth.usa.its.pdf")
varImpPlot(rf.study_group)
dev.off()
# Nigeria
asv_tab2 <- asv_tab[!grepl("L", row.names(asv_tab)),] #remove USA samples
asv_tab_norm <- sweep(asv_tab2, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata2 <- metadata[!grepl("L", row.names(metadata)),] #remove USA samples
metadata2$Tooth_Classification <- factor(metadata2$Tooth_Classification)
asv_tab_var$var <- metadata2[rownames(asv_tab_var), "Tooth_Classification"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.tooth.nigeria.its.pdf")
varImpPlot(rf.study_group)
dev.off()


#within location mouth health
# USA
asv_tab2 <- asv_tab[!grepl("DM00", row.names(asv_tab)),] #remove nigeria samples
asv_tab_norm <- sweep(asv_tab2, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata2 <- metadata[!grepl("DM00", row.names(metadata)),] #remove nigeria samples
metadata2$mouth_health <- factor(metadata2$mouth_health)
asv_tab_var$var <- metadata2[rownames(asv_tab_var), "mouth_health"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.mouth.usa.its.pdf")
varImpPlot(rf.study_group)
dev.off()
# Nigeria
asv_tab2 <- asv_tab[!grepl("L", row.names(asv_tab)),] #remove USA samples
asv_tab_norm <- sweep(asv_tab2, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata2 <- metadata[!grepl("L", row.names(metadata)),] #remove USA samples
metadata2$mouth_health <- factor(metadata2$mouth_health)
asv_tab_var$var <- metadata2[rownames(asv_tab_var), "mouth_health"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.mouth.nigeria.its.pdf")
varImpPlot(rf.study_group)
dev.off()

#within location tooth health
# USA
asv_tab2 <- asv_tab[!grepl("DM00", row.names(asv_tab)),] #remove nigeria samples
asv_tab_norm <- sweep(asv_tab2, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata2 <- metadata[!grepl("DM00", row.names(metadata)),] #remove nigeria samples
metadata2$tooth_health <- factor(metadata2$tooth_health)
asv_tab_var$var <- metadata2[rownames(asv_tab_var), "tooth_health"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.tooth_health.usa.its.pdf")
varImpPlot(rf.study_group)
dev.off()
# Nigeria
asv_tab2 <- asv_tab[!grepl("L", row.names(asv_tab)),] #remove USA samples
asv_tab_norm <- sweep(asv_tab2, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata2 <- metadata[!grepl("L", row.names(metadata)),] #remove USA samples
metadata2$tooth_health <- factor(metadata2$tooth_health)
asv_tab_var$var <- metadata2[rownames(asv_tab_var), "tooth_health"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.tooth_health.nigeria.its.pdf")
varImpPlot(rf.study_group)
dev.off()

#tooth health
asv_tab_norm <- sweep(asv_tab, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata$tooth_health <- factor(metadata$tooth_health)
asv_tab_var$var <- metadata[rownames(asv_tab_var), "tooth_health"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.tooth_health.its.pdf")
varImpPlot(rf.study_group)
dev.off()

#mouth health
asv_tab_norm <- sweep(asv_tab, 2, colSums(asv_tab2), '/')*100
asv_tab_scale <- scale(asv_tab_norm, center=T, scale=T)
asv_tab_var <- data.frame(asv_tab_scale)
asv_tab_var <- asv_tab_var[, colSums(is.na(asv_tab_var)) == 0]#remove NAs
metadata$mouth_health <- factor(metadata$mouth_health)
asv_tab_var$var <- metadata[rownames(asv_tab_var), "mouth_health"]
rf.study_group <- randomForest(x=asv_tab_var[,1:(ncol(asv_tab_var)-1)], y=asv_tab_var$var, ntree=10000, importance=T, proximity=T)
rf.study_group
pdf("./rf.mouth.its.pdf")
varImpPlot(rf.study_group)
dev.off()