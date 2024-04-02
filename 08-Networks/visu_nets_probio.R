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
  vertex.color=vcol, vertex.label=V(net.usa)$genus, vertex.size=vs)
dev.off()

pdf("test.pdf")
plot(net.usa ,mark.groups = split(V(net.usa)$name,V(net.usa)$genus))
dev.off()


#wc <- cluster_leading_eigen(net.usa, weights =NA)
# test.net <- delete.edges(net.usa, E(net.usa)[E(net.usa)$weight >0.3])
wc <- cluster_fast_greedy(net.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
#df %>%  filter(names == "rASV2" | names == "ASV4" | names == "rASV400" |names ==  "rASV864" |names ==  "rASV117"| names == "rASV105" | names =="rASV1327")
#df %>%  filter(x ==5) %>%left_join(genus, by=c('names'='ASV'))
#df %>%left_join(genus, by=c('names'='ASV')) %>% filter(V8 == "Streptococcus_parasanguinis"| V8 == "Limosilactobacillus_fermentum" | V8 == "Streptococcus_mutans" | 
#	V8 == "Scardovia_wiggsiae"| V8 == "Rhodotorula_mucilaginosa" | V8 == "Debaryomyces_prosopidis" | V8 =="Leptotrichia_sp._oral_taxon_212")

df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | 
  V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum")

df2 <- df %>%  filter(x ==8) %>%left_join(genus, by=c('names'='ASV'))
cluster_interest <- as.vector(df2$names)
my_fun <- function(i) {
 as.numeric(V(net.usa)[i])
}
usa_cluster <- sapply(cluster_interest, my_fun)

cluster.eight.usa <- induced_subgraph(net.usa, usa_cluster)
wc <- cluster_fast_greedy(cluster.eight.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | 
  V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum")

#nwc <- cluster_leading_eigen(net.nigeria, weights =NA)
#test.net <- delete.edges(net.nigeria, E(net.nigeria)[E(net.nigeria)$weight >0.3])
nwc <- cluster_fast_greedy(net.nigeria, weights =NA)
df <- as.data.frame(membership(nwc))
df$names <- rownames(df)
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | 
  V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum")

df2 <- df %>%  filter(x ==6) %>%left_join(genus, by=c('names'='ASV'))
cluster_interest <- as.vector(df2$names)
my_fun <- function(i) {
 as.numeric(V(net.nigeria)[i])
}
nigeria_cluster <- sapply(cluster_interest, my_fun)

cluster.six.nigeria <- induced_subgraph(net.nigeria, nigeria_cluster)
wc <- cluster_fast_greedy(cluster.six.nigeria, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | 
  V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | 
  V8 == "Corynebacterium_durum" | V8 =="Streptococcus_mutans")
df %>%  filter(x ==4) %>%left_join(genus, by=c('names'='ASV'))
