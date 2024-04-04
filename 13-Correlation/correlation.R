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
library(Hmisc)
library(tidyverse)
library(corrplot)
library(tabletools)
#set seed
set.seed(12349)
setwd("~/usa_nigeria/correlation")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

#transform data
its.dat.clr <- microbiome::transform(tax_glom(its.dat, "V8"), transform="clr", target="OTU")
rpoc.dat.clr <- microbiome::transform(tax_glom(rpoc.dat, "V8"), transform="clr", target="OTU")
#combine data
ps.dat.clr <- merge_phyloseq(rpoc.dat.clr, its.dat.clr)
#get frequency table
ASV_freq <- as(otu_table(ps.dat.clr), "matrix")
trans_ASV_freq <- t(ASV_freq)
trans_ASV_freq_df1 <- as.data.frame(trans_ASV_freq)
trans_ASV_freq_df2 <- rownames_to_column(trans_ASV_freq_df1, var = "ASV")
#getting taxa species
taxa = as(tax_table(ps.dat.clr)), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
#renmaing all to species level
renamed_ASV_freq1 <- left_join(trans_ASV_freq_df2, orderdf, by=c('ASV'='ASV'))
colnames(renamed_ASV_freq1)[113] <- "Species"
renamed_ASV_freq2 <- renamed_ASV_freq1[, c(113,2:112)]
renamed_ASV_freq4 <- as.data.frame(renamed_ASV_freq2)
rownames(renamed_ASV_freq4) <- renamed_ASV_freq4[,1]
renamed_ASV_freq4[,1] <- NULL
renamed_ASV_freq5 <- t(as.matrix(renamed_ASV_freq4))


#running correlation
its_rpoc_corr <- rcorr(renamed_ASV_freq5,type=c("spearman"))
#adjust p value
cor_enviro_adjust <- rcorr_padjust(its_rpoc_corr, method = "BH")
diag(cor_enviro_adjust$P) <- 0

pdf("cross_corr.pdf", width=100, height=100)
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off()

library(MASS)
library(vcd)
usa_meta <- filter(metadata, !Race == "NAN")
test <- assocstats(xtabs(~usa_meta$Ethnicity + usa_meta$Tooth_Classification))
p.adjust(test, method = "fdr", n = length(p))

assocstats(xtabs(~metadata$Geog_loc + metadata$Tooth_Classification))