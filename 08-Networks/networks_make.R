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
links <- subset(adjacencyData.usa, weight > 0) #only positive
#replace asv with spiecieas
nodes <- select(adjacencyData.usa, to) %>% unique()
nodes2 <- left_join(nodes, genus, by=c('to'='ASV'))
colnames(nodes2)[2] <- "genus"
nodes3 <- left_join(nodes2, domain, by=c('to'='ASV'))
colnames(nodes3)[3] <- "domain"
colnames(nodes3)[1] <- "ASV"

#convert back to igraph
net.usa <- graph_from_data_frame(d=links, vertices=nodes3, directed=T)

import_taxa = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa","Stereum_hirsutum",
	"Streptococcus_mutans", "Treponema_denticola", 
 "Corynebacterium_durum","Corynebacterium_segmentosum", "Corynebacterium_mustelae", "Corynebacterium_matruchotii")
#edges color
ecol <- rep("gray80", ecount(net.usa))
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Debaryomyces_prosopidis"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#332288"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Candida_tropicalis"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#117733"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Candida_albicans"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#CC6677"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Candida_unknown"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#88CCEE"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Rhodotorula_mucilaginosa"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#999933"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Stereum_hirsutum"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#882255"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Streptococcus_mutans"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#44AA99"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Treponema_denticola"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#DDCC77"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Corynebacterium_durum"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Corynebacterium_segmentosum"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Corynebacterium_mustelae"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Corynebacterium_matruchotii"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"

inc.edges <- incident(net.usa,  V(net.usa)[genus=="Debaryomyces_prosopidis"], mode="all")
ecol[inc.edges] <- "#332288"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Candida_tropicalis"], mode="all")
ecol[inc.edges] <- "#117733"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Candida_albicans"], mode="all")
ecol[inc.edges] <- "#CC6677"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Candida_unknown"], mode="all")
ecol[inc.edges] <- "#88CCEE"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Rhodotorula_mucilaginosa"], mode="all")
ecol[inc.edges] <- "#999933"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Stereum_hirsutum"], mode="all")
ecol[inc.edges] <- "#882255"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Streptococcus_mutans"], mode="all")
ecol[inc.edges] <- "#44AA99"
#inc.edges <- incident(net.usa,  V(net.usa)[genus=="Treponema_denticola"], mode="all")
#ecol[inc.edges] <- "#DDCC77"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Corynebacterium_durum"], mode="all")
ecol[inc.edges] <- "#AA4499"
#inc.edges <- incident(net.usa,  V(net.usa)[genus=="Corynebacterium_segmentosum"], mode="all")
#ecol[inc.edges] <- "#AA4499"
#inc.edges <- incident(net.usa,  V(net.usa)[genus=="Corynebacterium_mustelae"], mode="all")
#ecol[inc.edges] <- "#AA4499"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Corynebacterium_matruchotii"], mode="all")
ecol[inc.edges] <- "#AA4499"

#edge width
import_taxa2 = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa","Stereum_hirsutum",
	"Streptococcus_mutans","Corynebacterium_durum", "Corynebacterium_matruchotii")
ew <- rep(2, ecount(net.usa))
for(i in import_taxa) {
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus==i], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ew[inc.edges] <- 10
}
for(i in import_taxa2) {
inc.edges <- incident(net.usa,  V(net.usa)[genus==i], mode="all")
ew[inc.edges] <- 30
}
#vertix color
vcol <- adjustcolor(c("gray25", "gray60")[1+(V(net.usa)$domain=="Bacteria")], alpha=.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Debaryomyces_prosopidis"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#332288"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Candida_tropicalis"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#117733"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Candida_albicans"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#CC6677"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Candida_unknown"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#88CCEE"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Rhodotorula_mucilaginosa"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#999933"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Stereum_hirsutum"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#882255"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Streptococcus_mutans"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#44AA99"), alpha =.6)
#neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Treponema_denticola"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#DDCC77"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Corynebacterium_durum"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
#neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Corynebacterium_segmentosum"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
#neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Corynebacterium_mustelae"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Corynebacterium_matruchotii"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)

vcol[V(net.usa)$genus=="Debaryomyces_prosopidis"] <-  adjustcolor(c("#332288"), alpha =.9)
vcol[V(net.usa)$genus=="Candida_tropicalis"] <-  adjustcolor(c("#117733"), alpha =.9)
vcol[V(net.usa)$genus=="Candida_albicans"] <-  adjustcolor(c("#CC6677"), alpha =.9)
vcol[V(net.usa)$genus=="Candida_unknown"] <-  adjustcolor(c("#88CCEE"), alpha =.9)
vcol[V(net.usa)$genus=="Rhodotorula_mucilaginosa"] <-  adjustcolor(c("#999933"), alpha =.9)
vcol[V(net.usa)$genus=="Stereum_hirsutum"] <-  adjustcolor(c("#882255"), alpha =.9)
vcol[V(net.usa)$genus=="Streptococcus_mutans"] <-  adjustcolor(c("#44AA99"), alpha =.9)
vcol[V(net.usa)$genus=="Treponema_denticola"] <-  adjustcolor(c("#DDCC77"), alpha =.9)
vcol[V(net.usa)$genus=="Corynebacterium_durum"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.usa)$genus=="Corynebacterium_segmentosum"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.usa)$genus=="Corynebacterium_mustelae"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.usa)$genus=="Corynebacterium_matruchotii"] <-  adjustcolor(c("#AA4499"), alpha =.9)
#vertix border color
vbod <- adjustcolor(c("gray15", "gray70")[1+(V(net.usa)$domain=="Bacteria")], alpha=.8)
	
#vertix size
vs <- rep(3, vcount(net.usa))

for(i in import_taxa) {
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus==i], mode="all", order = 1))
vs[neigh.nodes] <- 4
}

for(i in import_taxa2) {
vs[V(net.usa)$genus==i] <- 10
}
#names
font <- rep(.001, vcount(net.usa))

for(i in import_taxa) {
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus==i], mode="all", order = 1))
font[neigh.nodes] <- 2.5
}

for(i in import_taxa2) {
font[V(net.usa)$genus==i] <- 5
}
#font color
fontc <- rep(adjustcolor("black", alpha = 0.0), vcount(net.usa))
for(i in import_taxa2) {
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus==i], mode="all")
fontc[neigh.nodes] <- "black"
}
for(i in import_taxa2) {
fontc[V(net.usa)$genus==i] <- "black"
}
#font type
fontt <- rep(1, vcount(net.usa))

for(i in import_taxa2) {
fontt[V(net.usa)$genus==i] <- 2
}
#remove fungi vertixs
net.usa.bact <- subgraph.edges(graph=net.usa, eids=which(V(net.usa)$domain=="Fungi"), delete.vertices = TRUE)







#nigeria
getStability(nigeria.cross)
#nigeria.cross
sebeta.nigeria.cross <- symBeta(getOptBeta(nigeria.cross), mode='maxabs')
taxnames <- c(taxa_names(rpoc.nigeria), taxa_names(its.nigeria))
colnames(sebeta.nigeria.cross) <- rownames(sebeta.nigeria.cross) <- taxnames
weighted.adj.nigeria.cross <- sebeta.nigeria.cross*getRefit(nigeria.cross)
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
#nigeria
spiec_df_nigeria <- as.matrix(weighted.adj.nigeria.cross)
adjacencyData.nigeria <- as.data.frame(melt(spiec_df_nigeria))
colnames(adjacencyData.nigeria)[1] <- "to"
colnames(adjacencyData.nigeria)[2] <- "from"
colnames(adjacencyData.nigeria)[3] <- "weight"
#making links
links <- subset(adjacencyData.nigeria, weight > 0) #only positive
#replace asv with spiecieas
nodes <- select(adjacencyData.nigeria, to) %>% unique()
nodes2 <- left_join(nodes, genus, by=c('to'='ASV'))
colnames(nodes2)[2] <- "genus"
nodes3 <- left_join(nodes2, domain, by=c('to'='ASV'))
colnames(nodes3)[3] <- "domain"
colnames(nodes3)[1] <- "ASV"

#convert back to igraph
net.nigeria <- graph_from_data_frame(d=links, vertices=nodes3, directed=T)

import_taxa = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa","Stereum_hirsutum",
	"Streptococcus_mutans", "Treponema_denticola", 
 "Corynebacterium_durum","Corynebacterium_segmentosum", "Corynebacterium_mustelae", "Corynebacterium_matruchotii")
#edges color
ecol <- rep("gray80", ecount(net.nigeria))
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Debaryomyces_prosopidis"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#332288"
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Candida_tropicalis"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#117733"
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Candida_albicans"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#CC6677"
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Candida_unknown"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#88CCEE"
#neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Rhodotorula_mucilaginosa"], mode="all", order = 1))
#inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
#ecol[inc.edges] <- "#999933"
#neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Stereum_hirsutum"], mode="all", order = 1))
#inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
#ecol[inc.edges] <- "#882255"
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Streptococcus_mutans"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#44AA99"
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Treponema_denticola"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#DDCC77"
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_durum"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
#neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_segmentosum"], mode="all", order = 1))
#inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
#ecol[inc.edges] <- "#AA4499"
#neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_mustelae"], mode="all", order = 1))
#inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
#ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_matruchotii"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"

inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Debaryomyces_prosopidis"], mode="all")
ecol[inc.edges] <- "#332288"
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Candida_tropicalis"], mode="all")
ecol[inc.edges] <- "#117733"
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Candida_albicans"], mode="all")
ecol[inc.edges] <- "#CC6677"
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Candida_unknown"], mode="all")
ecol[inc.edges] <- "#88CCEE"
#inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Rhodotorula_mucilaginosa"], mode="all")
#ecol[inc.edges] <- "#999933"
#inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Stereum_hirsutum"], mode="all")
#ecol[inc.edges] <- "#882255"
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Streptococcus_mutans"], mode="all")
ecol[inc.edges] <- "#44AA99"
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Treponema_denticola"], mode="all")
ecol[inc.edges] <- "#DDCC77"
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Corynebacterium_durum"], mode="all")
ecol[inc.edges] <- "#AA4499"
#inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Corynebacterium_segmentosum"], mode="all")
#ecol[inc.edges] <- "#AA4499"
#inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Corynebacterium_mustelae"], mode="all")
#ecol[inc.edges] <- "#AA4499"
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus=="Corynebacterium_matruchotii"], mode="all")
ecol[inc.edges] <- "#AA4499"

#edge width
ew <- rep(2, ecount(net.nigeria))
import_taxa3 = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Streptococcus_mutans", 
	"Treponema_denticola", "Corynebacterium_durum", "Corynebacterium_matruchotii")
for(i in import_taxa3) {
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus==i], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.nigeria, neigh.nodes, mode = "all"))
ew[inc.edges] <- 10
}
for(i in import_taxa3) {
inc.edges <- incident(net.nigeria,  V(net.nigeria)[genus==i], mode="all")
ew[inc.edges] <- 30
}
#vertix color
vcol <- adjustcolor(c("gray25", "gray60")[1+(V(net.nigeria)$domain=="Bacteria")], alpha=.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Debaryomyces_prosopidis"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#332288"), alpha =.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Candida_tropicalis"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#117733"), alpha =.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Candida_albicans"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#CC6677"), alpha =.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Candida_unknown"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#88CCEE"), alpha =.6)
#neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Rhodotorula_mucilaginosa"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#999933"), alpha =.6)
#neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Stereum_hirsutum"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#882255"), alpha =.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Streptococcus_mutans"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#44AA99"), alpha =.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Treponema_denticola"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#DDCC77"), alpha =.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_durum"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
#neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_segmentosum"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
#neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_mustelae"], mode="all")
#vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus=="Corynebacterium_matruchotii"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)

vcol[V(net.nigeria)$genus=="Debaryomyces_prosopidis"] <-  adjustcolor(c("#332288"), alpha =.9)
vcol[V(net.nigeria)$genus=="Candida_tropicalis"] <-  adjustcolor(c("#117733"), alpha =.7)
vcol[V(net.nigeria)$genus=="Candida_albicans"] <-  adjustcolor(c("#CC6677"), alpha =.9)
vcol[V(net.nigeria)$genus=="Candida_unknown"] <-  adjustcolor(c("#88CCEE"), alpha =.9)
vcol[V(net.nigeria)$genus=="Rhodotorula_mucilaginosa"] <-  adjustcolor(c("#999933"), alpha =.9)
vcol[V(net.nigeria)$genus=="Stereum_hirsutum"] <-  adjustcolor(c("#882255"), alpha =.9)
vcol[V(net.nigeria)$genus=="Streptococcus_mutans"] <-  adjustcolor(c("#44AA99"), alpha =.9)
vcol[V(net.nigeria)$genus=="Treponema_denticola"] <-  adjustcolor(c("#DDCC77"), alpha =.9)
vcol[V(net.nigeria)$genus=="Corynebacterium_durum"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.nigeria)$genus=="Corynebacterium_segmentosum"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.nigeria)$genus=="Corynebacterium_mustelae"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.nigeria)$genus=="Corynebacterium_matruchotii"] <-  adjustcolor(c("#AA4499"), alpha =.9)
#vertix border color
vbod <- adjustcolor(c("gray15", "gray70")[1+(V(net.nigeria)$domain=="Bacteria")], alpha=.8)
#vertix size
vs <- rep(3, vcount(net.nigeria))

for(i in import_taxa3) {
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus==i], mode="all", order = 1))
vs[neigh.nodes] <- 4
}

for(i in import_taxa3) {
vs[V(net.nigeria)$genus==i] <- 10
}
#names
font <- rep(.001, vcount(net.nigeria))

for(i in import_taxa3) {
neigh.nodes <- unlist(ego(net.nigeria, V(net.nigeria)[genus==i], mode="all", order = 1))
font[neigh.nodes] <- 2.5
}

for(i in import_taxa3) {
font[V(net.nigeria)$genus==i] <- 5
}
#font color
fontc <- rep(adjustcolor("black", alpha = 0.0), vcount(net.nigeria))
for(i in import_taxa3) {
neigh.nodes <- neighbors(net.nigeria, V(net.nigeria)[genus==i], mode="all")
fontc[neigh.nodes] <- "black"
}
for(i in import_taxa3) {
fontc[V(net.nigeria)$genus==i] <- "black"
}
#font type
fontt <- rep(1, vcount(net.nigeria))

for(i in import_taxa3) {
fontt[V(net.nigeria)$genus==i] <- 2
}

net.nigeria.bact <- subgraph.edges(graph=net.nigeria, eids=which(V(net.nigeria)$domain=="Fungi"), delete.vertices = TRUE)


save.image("networks.RData")
