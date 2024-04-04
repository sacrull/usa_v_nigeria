#devtools::install_github("JLSteenwyk/ggpubfigs")
library(ggpubfigs)
library(SpiecEasi)
library(igraph)
library(viridis)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(Matrix)
setwd("/home/suzanne/usa_nigeria/networks")
load("/home/suzanne/usa_nigeria/phyloseq_obj/ps.RData")
set.seed=123456
# USA rpoc
rpoc.usa <- subset_samples(rpoc.dat, Geog_loc == "USA")
glom.rpoc.usa <- tax_glom(rpoc.usa, "V8")
# filter low abundance taxa to simplify network building step (seen at least 3 times in 10% of samples)
rpoc.usa <- filter_taxa(glom.rpoc.usa, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
# USA its
its.usa <- subset_samples(its.dat, Geog_loc == "USA")
glom.its.usa <- tax_glom(its.usa, "V8")
# filter low abundance taxa to simplify network building step (seen at least 3 times in 10% of samples)
its.usa<- filter_taxa(glom.its.usa, function(x) sum(x > 3) > (0.05*length(x)), TRUE)

#Niger rpoc
rpoc.nigeria <- subset_samples(rpoc.dat, Geog_loc == "Nigeria")
glom.rpoc.nigeria <- tax_glom(rpoc.nigeria, "V8")
# filter low abundance taxa to simplify network building step (seen at least 3 times in 10% of samples)
rpoc.nigeria <- filter_taxa(glom.rpoc.nigeria, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
# Niger ITS
its.nigeria <- subset_samples(its.dat, Geog_loc == "Nigeria")
glom.its.nigeria <- tax_glom(its.nigeria, "V8")
# filter low abundance taxa to simplify network building step (seen at least 3 times in 10% of samples)
its.nigeria<- filter_taxa(glom.its.nigeria, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
#network generation
pargs <- list(rep.num = 100, ncores=7, seed=123456)


usa.cross <- spiec.easi(list(rpoc.usa, its.usa), method='mb', nlambda=75,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(usa.cross)

nigeria.cross <- spiec.easi(list(rpoc.nigeria, its.nigeria), method='mb', nlambda=75,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(nigeria.cross)

save.image("cross_nets.RData")
load("cross_nets.RData")
#usa.cross
sebeta.usa.cross <- symBeta(getOptBeta(usa.cross), mode='maxabs')
taxnames <- c(taxa_names(rpoc.usa), taxa_names(its.usa))
colnames(sebeta.usa.cross) <- rownames(sebeta.usa.cross) <- taxnames
weighted.adj.usa.cross <- sebeta.usa.cross*getRefit(usa.cross)
library(circlize)
library(tidyverse)
library(reshape2)
#getting just genus name
taxa = as(tax_table(ps.data), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
genus <- orderdf %>% 
  rownames_to_column(var = "ASV")
#getting phyla names
taxa = as(tax_table(ps.data), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V2)
domain <- orderdf %>% 
  rownames_to_column(var = "ASV")
#get orders
list_all<-as.list(orderdf)
ordered_list <- unique(sort(list_all$V8))
ordered_list2 <- c(ordered_list)

#making data frame from igraph
#usa
spiec_df_usa <- as.matrix(weighted.adj.usa.cross)
adjacencyData.usa <- as.data.frame(melt(spiec_df_usa))
colnames(adjacencyData.usa)[1] <- "to"
colnames(adjacencyData.usa)[2] <- "from"
colnames(adjacencyData.usa)[3] <- "weight"
#making links
links <- subset(adjacencyData.usa, weight < 0) #only negative
#replace asv with spiecieas
nodes <- select(adjacencyData.usa, to) %>% unique()
nodes2 <- left_join(nodes, genus, by=c('to'='ASV'))
colnames(nodes2)[2] <- "genus"
nodes3 <- left_join(nodes2, domain, by=c('to'='ASV'))
colnames(nodes3)[3] <- "domain"
colnames(nodes3)[1] <- "ASV"

#convert back to igraph
net <- graph_from_data_frame(d=links, vertices=nodes3, directed=T)
friendly_pal("muted_nine", 9, type = "discrete")

import_taxa = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa","Stereum_hirsutum",
	"Streptococcus_mutans", "Treponema_denticola", 
 "Corynebacterium_durum","Corynebacterium_segmentosum", "Corynebacterium_mustelae", "Corynebacterium_matruchotii")
#edges color
ecol <- rep("gray80", ecount(net))
neigh.nodes <- unlist(ego(net, V(net)[genus=="Debaryomyces_prosopidis"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#332288"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Candida_tropicalis"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#117733"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Candida_albicans"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#CC6677"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Candida_unknown"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#88CCEE"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Rhodotorula_mucilaginosa"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#999933"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Stereum_hirsutum"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#882255"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Streptococcus_mutans"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#44AA99"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Treponema_denticola"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#DDCC77"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Corynebacterium_durum"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Corynebacterium_segmentosum"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Corynebacterium_mustelae"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net, V(net)[genus=="Corynebacterium_matruchotii"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"

inc.edges <- incident(net,  V(net)[genus=="Debaryomyces_prosopidis"], mode="all")
ecol[inc.edges] <- "#332288"
inc.edges <- incident(net,  V(net)[genus=="Candida_tropicalis"], mode="all")
ecol[inc.edges] <- "#117733"
inc.edges <- incident(net,  V(net)[genus=="Candida_albicans"], mode="all")
ecol[inc.edges] <- "#CC6677"
inc.edges <- incident(net,  V(net)[genus=="Candida_unknown"], mode="all")
ecol[inc.edges] <- "#88CCEE"
inc.edges <- incident(net,  V(net)[genus=="Rhodotorula_mucilaginosa"], mode="all")
ecol[inc.edges] <- "#999933"
inc.edges <- incident(net,  V(net)[genus=="Stereum_hirsutum"], mode="all")
ecol[inc.edges] <- "#882255"
inc.edges <- incident(net,  V(net)[genus=="Streptococcus_mutans"], mode="all")
ecol[inc.edges] <- "#44AA99"
#inc.edges <- incident(net,  V(net)[genus=="Treponema_denticola"], mode="all")
#ecol[inc.edges] <- "#DDCC77"
inc.edges <- incident(net,  V(net)[genus=="Corynebacterium_durum"], mode="all")
ecol[inc.edges] <- "#AA4499"
#inc.edges <- incident(net,  V(net)[genus=="Corynebacterium_segmentosum"], mode="all")
#ecol[inc.edges] <- "#AA4499"
#inc.edges <- incident(net,  V(net)[genus=="Corynebacterium_mustelae"], mode="all")
#ecol[inc.edges] <- "#AA4499"
inc.edges <- incident(net,  V(net)[genus=="Corynebacterium_matruchotii"], mode="all")
ecol[inc.edges] <- "#AA4499"

#edge width
import_taxa2 = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa","Stereum_hirsutum",
	"Streptococcus_mutans","Corynebacterium_durum", "Corynebacterium_matruchotii")
ew <- rep(2, ecount(net))
for(i in import_taxa) {
neigh.nodes <- unlist(ego(net, V(net)[genus==i], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net, neigh.nodes, mode = "all"))
ew[inc.edges] <- 10
}
for(i in import_taxa2) {
inc.edges <- incident(net,  V(net)[genus==i], mode="all")
ew[inc.edges] <- 30
}
#vertix color
vcol <- adjustcolor(c("gray25", "gray60")[1+(V(net)$domain=="Bacteria")], alpha=.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Debaryomyces_prosopidis"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#332288"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Candida_tropicalis"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#117733"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Candida_albicans"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#CC6677"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Candida_unknown"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#88CCEE"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Rhodotorula_mucilaginosa"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#999933"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Stereum_hirsutum"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#882255"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Streptococcus_mutans"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#44AA99"), alpha =.6)
#neigh.nodes <- neighbors(net, V(net)[genus=="Treponema_denticola"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#DDCC77"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Corynebacterium_durum"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
#neigh.nodes <- neighbors(net, V(net)[genus=="Corynebacterium_segmentosum"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
#neigh.nodes <- neighbors(net, V(net)[genus=="Corynebacterium_mustelae"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
neigh.nodes <- neighbors(net, V(net)[genus=="Corynebacterium_matruchotii"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)

vcol[V(net)$genus=="Debaryomyces_prosopidis"] <-  adjustcolor(c("#332288"), alpha =.9)
vcol[V(net)$genus=="Candida_tropicalis"] <-  adjustcolor(c("#117733"), alpha =.9)
vcol[V(net)$genus=="Candida_albicans"] <-  adjustcolor(c("#CC6677"), alpha =.9)
vcol[V(net)$genus=="Candida_unknown"] <-  adjustcolor(c("#88CCEE"), alpha =.9)
vcol[V(net)$genus=="Rhodotorula_mucilaginosa"] <-  adjustcolor(c("#999933"), alpha =.9)
vcol[V(net)$genus=="Stereum_hirsutum"] <-  adjustcolor(c("#882255"), alpha =.9)
vcol[V(net)$genus=="Streptococcus_mutans"] <-  adjustcolor(c("#44AA99"), alpha =.9)
vcol[V(net)$genus=="Treponema_denticola"] <-  adjustcolor(c("#DDCC77"), alpha =.9)
vcol[V(net)$genus=="Corynebacterium_durum"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net)$genus=="Corynebacterium_segmentosum"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net)$genus=="Corynebacterium_mustelae"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net)$genus=="Corynebacterium_matruchotii"] <-  adjustcolor(c("#AA4499"), alpha =.9)
#vertix border color
vbod <- adjustcolor(c("gray15", "gray70")[1+(V(net)$domain=="Bacteria")], alpha=.8)
	
#vertix size
vs <- rep(3, vcount(net))

for(i in import_taxa) {
neigh.nodes <- unlist(ego(net, V(net)[genus==i], mode="all", order = 1))
vs[neigh.nodes] <- 4
}

for(i in import_taxa2) {
vs[V(net)$genus==i] <- 10
}
#names
font <- rep(.001, vcount(net))

for(i in import_taxa) {
neigh.nodes <- unlist(ego(net, V(net)[genus==i], mode="all", order = 1))
font[neigh.nodes] <- 2.5
}

for(i in import_taxa2) {
font[V(net)$genus==i] <- 5
}
#font color
fontc <- rep(adjustcolor("black", alpha = 0.0), vcount(net))
for(i in import_taxa2) {
neigh.nodes <- neighbors(net, V(net)[genus==i], mode="all")
fontc[neigh.nodes] <- "black"
}
for(i in import_taxa2) {
fontc[V(net)$genus==i] <- "black"
}
#font type
fontt <- rep(1, vcount(net))

for(i in import_taxa2) {
fontt[V(net)$genus==i] <- 2
}

#graph
pdf("usa_negative.pdf", height = 100, width = 100)
plot(net,vertex.color=vcol, edge.color=ecol, edge.width=ew, vertex.label=V(net)$genus, vertex.size=vs, vertex.frame.color =vbod,
 edge.arrow.mode=0, vertex.frame.color = "white", vertex.label.cex = font, vertex.label.color = fontc, vertex.label.font=fontt)
dev.off()

ceb <- cluster_edge_betweenness(net) 
pdf("test.pdf", height = 100, width = 100)
dendPlot(ceb, mode="hclust")
dev.off()



