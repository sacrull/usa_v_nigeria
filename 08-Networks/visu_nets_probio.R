library(brainGraph)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(scales)
library(phyloseq)
library(SpiecEasi)
library(ggpubr)
library(maxnodf)
library(RColorBrewer)
library(ggside)
library(data.table)
library(foreach)
setwd("/home/suzanne/usa_nigeria/networks/probio")
load("/home/suzanne/usa_nigeria/networks/probio/cross_net.RData")

#vertix colors
vcol <- adjustcolor(c("gray25", "#F4E6BF")[1+(V(net.usa)$domain=="Bacteria")], alpha=.6)
#vertix size
vs <- rep(3, vcount(net.usa))

wc <- cluster_leading_eigen(net.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
cluster <- as.vector(membership(wc))

df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"|
 V8 == "Streptococcus_oralis" | V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | 
 V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum")

cluster_colors <- ifelse(cluster == 4, "green", ifelse(cluster == 6, "purple", adjustcolor(c("#AA4499"), alpha =.3)))

length(wc)
adjustcolor(colors, alpha =.3)
colors <- c("gray", "gray", "gray", "purple", "gray","green", "gray", "gray", "blue", "gray", "gray", "gray", "gray")
pdf("test.pdf", width =100, height =100)
plot(net.usa ,mark.groups = wc , mark.shape = 1/2,
  mark.col = adjustcolor(colors, alpha =.3),
  vertex.color= membership(wc), vertex.label=V(net.usa)$genus, vertex.size=vs)
dev.off()

pdf("test.pdf")
plot(net.usa ,mark.groups = split(V(net.usa)$name,V(net.usa)$genus))
dev.off()




#cluster network
wc <- cluster_leading_eigen(net.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
#get clusters that have species of interest
import_clusts <- df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum" | V8 == "Lautropia_mirabilis" | V8 == "Haemophilus_parainfluenzae" | V8 == "Haemophilus_influenzae")
import_clusts2 <- unique(as.vector(import_clusts$x))
important_clusters <- function(i) {
	filter(df, x ==i)
}
import_clust3 <- lapply(import_clusts2, important_clusters)
#get ASVs that are in those clusters
import_clust4 <- as.data.frame(unlist(import_clust3))
colnames(import_clust4)[1] <- "names"
df2 <- as.data.frame(import_clust4[grep("ASV", import_clust4$names), ])
colnames(df2)[1] <- "names"
cluster_interest <- as.vector(df2$names)
my_fun <- function(i) {
 as.numeric(V(net.usa)[i])
}
usa_cluster <- sapply(cluster_interest, my_fun)
#subset graph with clusters
cluster.usa <- induced_subgraph(net.usa, usa_cluster)
#recluster network
wc <- cluster_leading_eigen(cluster.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Rhodotorula_mucilaginosa" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum" | V8 == "Lautropia_mirabilis" | V8 == "Haemophilus_parainfluenzae" | V8 == "Haemophilus_influenzae" | V8 == "Candida_albicans" | V8 == "Streptococcus_mutans")
df %>%  filter(x ==7) %>%left_join(genus, by=c('names'='ASV'))
#cluster color
length(wc)
colors <- c("#BEBEBE00", adjustcolor("#DD77F5", alpha =.5), adjustcolor("#DD7753", alpha =.5), adjustcolor("#3359C2", alpha =.5), "#BEBEBE00", adjustcolor("#68B160", alpha =.5), "#BEBEBE00", adjustcolor("#78BEFC", alpha =.5))
border <- c("#BEBEBE00", "#DD77F5", "#DD7753", "#3359C2", "#BEBEBE00", "#68B160", "#BEBEBE00", "#78BEFC")
#vertix colors
# vcol <- adjustcolor(c("gray25", "#F4E6BF")[1+(V(cluster.usa)$domain=="Bacteria")], alpha=.6)
v <- membership(wc)
vcol = ifelse(v == 2, "#DD77F5" , ifelse( v== 3, "#DD7753", ifelse(v == 4, "#3359C2", ifelse( v==6, "#68B160", ifelse( v ==8, "#78BEFC", "#F4E6BF")))))
names(vcol) <- NULL
#vertix size
import_taxa2 <- c("Streptococcus_sanguinis", "Streptococcus_oralis", "Leptotrichia_sp._oral_taxon_212", "Corynebacterium_durum", "Leptotrichia_sp._oral_taxon_215", "Lautropia_mirabilis", "Haemophilus_influenzae","Haemophilus_parainfluenzae")
vs <- rep(3, vcount(cluster.usa))
for(i in import_taxa2) {
vs[V(cluster.usa)$genus==i] <- 7
}
#vertix shape
vss <- c("square", "circle")[1+(V(cluster.usa)$domain=="Bacteria")]
#make network
modularity(wc)
pdf("usa_probionet.pdf", width =50, height =50)
plot(cluster.usa, layout= layout_with_gem, mark.groups = wc , mark.shape = 1/2,mark.col =colors, mark.border = border, 
	vertex.color = vcol,vertex.label=V(cluster.usa)$genus, vertex.size=vs, vertex.shape =vss)
dev.off()

pdf("test.pdf")
plot(cluster.usa,layout= layout_with_graphopt,mark.groups = wc, vertex.color = vcol,
 vertex.label=V(cluster.usa)$genus, vertex.size=vs, vertex.shape =vss)
dev.off()








#nigeria
wc <- cluster_leading_eigen(net.nigeria, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
#get clusters that have species of interest
import_clusts <- df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum" | V8 == "Lautropia_mirabilis" | V8 == "Haemophilus_parainfluenzae" | V8 == "Haemophilus_influenzae")
import_clusts2 <- unique(as.vector(import_clusts$x))
important_clusters <- function(i) {
	filter(df, x ==i)
}
import_clust3 <- lapply(import_clusts2, important_clusters)
#get ASVs that are in those clusters
import_clust4 <- as.data.frame(unlist(import_clust3))
colnames(import_clust4)[1] <- "names"
df2 <- as.data.frame(import_clust4[grep("ASV", import_clust4$names), ])
colnames(df2)[1] <- "names"
cluster_interest <- as.vector(df2$names)
my_fun <- function(i) {
 as.numeric(V(net.nigeria)[i])
}
nigiera_cluster <- sapply(cluster_interest, my_fun)
#subset graph with clusters
cluster.nigeria <- induced_subgraph(net.nigeria, nigiera_cluster)
wc <- cluster_leading_eigen(cluster.nigeria, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum" | V8 == "Candida_albicans" | V8 == "Streptococcus_mutans"| V8 == "Lautropia_mirabilis" | V8 == "Haemophilus_parainfluenzae" | V8 == "Haemophilus_influenzae")

df %>%  filter(x ==1) %>%left_join(genus, by=c('names'='ASV'))

#cluster color
length(wc)
colors <- c("#BEBEBE00", adjustcolor( "#DD77F5", alpha =.5),  adjustcolor("#68B160", alpha =.5), adjustcolor("#3359C2", alpha =.5), "#BEBEBE00", "#BEBEBE00", "#BEBEBE00", adjustcolor("#DD7753", alpha =.5))
border <- c("#BEBEBE00","#DD77F5", "#68B160", "#3359C2", "#BEBEBE00", "#BEBEBE00", "#BEBEBE00", "#DD7753")
#vertix colors
# vcol <- adjustcolor(c("gray25", "#F4E6BF")[1+(V(cluster.usa)$domain=="Bacteria")], alpha=.6)
v <- membership(wc)
vcol = ifelse(v == 2, "#DD77F5" , ifelse( v== 3, "#68B160", ifelse( v== 4, "#3359C2", ifelse( v== 8, "#DD7753", "#F4E6BF"))))
names(vcol) <- NULL
#vertix size
import_taxa2 <- c("Streptococcus_sanguinis", "Streptococcus_oralis", "Leptotrichia_sp._oral_taxon_212", "Corynebacterium_durum", "Leptotrichia_sp._oral_taxon_215", "Streptococcus_parasanguinis")
vs <- rep(3, vcount(cluster.nigeria))
for(i in import_taxa2) {
vs[V(cluster.nigeria)$genus==i] <- 7
}
#vertix shape
vss <- c("square", "circle")[1+(V(cluster.nigeria)$domain=="Bacteria")]
modularity(wc)
pdf("nigeria_probionet.pdf", width =50, height =50)
plot(cluster.nigeria, layout=layout_with_graphopt, mark.groups = wc , mark.shape = 1/2,mark.col =colors, mark.border = border, 
	vertex.color = vcol,vertex.label=V(cluster.nigeria)$genus, vertex.size=vs, vertex.shape =vss)
dev.off()