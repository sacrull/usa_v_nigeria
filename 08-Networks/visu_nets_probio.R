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
vcol <- adjustcolor(c("gray25", "gray60")[1+(V(net.usa)$domain=="Bacteria")], alpha=.6)
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





wc <- cluster_leading_eigen(net.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
#get clusters that have species of interest
import_clusts <- df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum" )
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
wc <- cluster_leading_eigen(cluster.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum" | V8 == "Candida_albicans" | V8 == "Streptococcus_mutans")


df %>%  filter(x ==7) %>%left_join(genus, by=c('names'='ASV'))
#cluster color
length(wc)
colors <- c(adjustcolor("#406E80", alpha =.5), "#BEBEBE00", adjustcolor("#DD7753", alpha =.5), "#BEBEBE00", "#BEBEBE00","#BEBEBE00", adjustcolor("#3359C2", alpha =.5), "#BEBEBE00", adjustcolor("#68B160", alpha =.5))
border <- c("#406E80", "gray60", "#DD7753", "gray60", "gray60","gray60", "#3359C2", "gray", "#68B160")
#vertix colors
# vcol <- adjustcolor(c("gray25", "gray60")[1+(V(cluster.usa)$domain=="Bacteria")], alpha=.6)
v <- membership(wc)
vcol = ifelse(v == 1, "#406E80" , ifelse( v== 3, "#DD7753", ifelse(v == 7, "#3359C2", ifelse( v==9, "#68B160", "gray60"))))
names(vcol) <- NULL
#vertix size
import_taxa2 <- c("Streptococcus_sanguinis", "Streptococcus_oralis", "Leptotrichia_sp._oral_taxon_212", "Corynebacterium_durum", "Leptotrichia_sp._oral_taxon_215", "Streptococcus_parasanguinis")
vs <- rep(3, vcount(cluster.usa))
for(i in import_taxa2) {
vs[V(cluster.usa)$genus==i] <- 7
}
#vertix shape
vss <- c("square", "circle")[1+(V(cluster.usa)$domain=="Bacteria")]

pdf("test.pdf", width =100, height =100)
plot(cluster.usa , mark.groups = wc , mark.shape = 1/2,mark.col =colors, mark.border = border, vertex.color = vcol,
 vertex.label=V(cluster.usa)$genus, vertex.size=vs, vertex.shape =vss)
dev.off()











#nigeria
wc <- cluster_leading_eigen(net.nigeria, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
#get clusters that have species of interest
import_clusts <- left_join(df, genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" )
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
usa_cluster <- sapply(cluster_interest, my_fun)
#subset graph with clusters
cluster.nigeria <- induced_subgraph(net.nigeria, usa_cluster)
wc <- cluster_leading_eigen(cluster.nigeria, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum" | V8 == "Candida_albicans" | V8 == "Streptococcus_mutans")

df %>%  filter(x ==1) %>%left_join(genus, by=c('names'='ASV'))


length(wc)
colors <- c(adjustcolor("purple", alpha =.3), "#BEBEBE00", "#BEBEBE00", "#BEBEBE00", "#BEBEBE00","#BEBEBE00", "#BEBEBE00", "#BEBEBE00", "#BEBEBE00", "#BEBEBE00", "#BEBEBE00", "#BEBEBE00", "#BEBEBE00")
adjustcolor(colors, alpha =.3)

vcol <- adjustcolor(c("gray25", "gray60")[1+(V(cluster.nigeria)$domain=="Bacteria")], alpha=.6)

pdf("test.pdf", width =100, height =100)
plot(cluster.nigeria , mark.groups = wc , mark.shape = 1,mark.col =colors,
  vertex.color=vcol, vertex.label=V(cluster.nigeria)$genus, vertex.size=3)
dev.off()