library(tidyverse)
library(reshape2)
setwd("/home/suzanne/usa_nigeria/networks/probio")
load("cross_nets.RData") 
#usa.cross
sebeta.usa.cross <- symBeta(getOptBeta(usa.cross), mode='maxabs')
taxnames <- c(taxa_names(rpoc.usa), taxa_names(its.usa))
colnames(sebeta.usa.cross) <- rownames(sebeta.usa.cross) <- taxnames
weighted.adj.usa.cross <- sebeta.usa.cross*getRefit(usa.cross)

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
net.usa <- adj2igraph(weighted.adj.usa.cross, rmEmptyNodes = TRUE, diag = FALSE,
  edge.attr = list(), vertex.attr = nodes3)
net.usa <- delete.edges(net.usa,which(E(net.usa)$weight<0))

#usa bact
sebeta.usa.cross <- symBeta(getOptBeta(usa.bact), mode='maxabs')
taxnames <- c(taxa_names(rpoc.usa))
colnames(sebeta.usa.cross) <- rownames(sebeta.usa.cross) <- taxnames
weighted.adj.usa.cross <- sebeta.usa.cross*getRefit(usa.bact)
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
net.usa.bact <- adj2igraph(weighted.adj.usa.cross, rmEmptyNodes = TRUE, diag = FALSE,
  edge.attr = list(), vertex.attr = nodes3)
net.usa.bact <- delete.edges(net.usa.bact,which(E(net.usa.bact)$weight<0))

#nigeria.cross
sebeta.nigeria.cross <- symBeta(getOptBeta(nigeria.cross), mode='maxabs')
taxnames <- c(taxa_names(rpoc.nigeria), taxa_names(its.nigeria))
colnames(sebeta.nigeria.cross) <- rownames(sebeta.nigeria.cross) <- taxnames
weighted.adj.nigeria.cross <- sebeta.nigeria.cross*getRefit(nigeria.cross)
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
net.nigeria <- adj2igraph(weighted.adj.nigeria.cross, rmEmptyNodes = TRUE, diag = FALSE,
  edge.attr = list(), vertex.attr = nodes3)
net.nigeria <- delete.edges(net.nigeria,which(E(net.nigeria)$weight<0))

#nigeria bact
sebeta.usa.bact <- symBeta(getOptBeta(nigeria.bact), mode='maxabs')
taxnames <- c(taxa_names(rpoc.nigeria))
colnames(sebeta.usa.bact) <- rownames(sebeta.usa.bact) <- taxnames
weighted.adj.nigeria.bact <- sebeta.usa.bact*getRefit(nigeria.bact)

#making data frame from igraph
spiec_df_usa <- as.matrix(weighted.adj.nigeria.bact)
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
net.nigeria.bact <- adj2igraph(weighted.adj.nigeria.bact, rmEmptyNodes = TRUE, diag = FALSE,
  edge.attr = list(), vertex.attr = nodes3)
net.nigeria.bact <- delete.edges(net.nigeria.bact,which(E(net.nigeria.bact)$weight<0))

save.image("networks_2.RData")