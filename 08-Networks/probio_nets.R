#devtools::install_github("JLSteenwyk/ggpubfigs")
library(ggpubfigs)
library(SpiecEasi)
library(igraph)
library(viridis)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(Matrix)
setwd("/home/suzanne/usa_nigeria/networks/probio")
load("/home/suzanne/usa_nigeria/phyloseq_obj/ps.RData")
set.seed=96789
# USA rpoc
rpoc.usa <- subset_samples(rpoc.dat, Geog_loc == "USA")
glom.rpoc.usa <- tax_glom(rpoc.usa, "V8")
# filter low abundance taxa to simplify net.usawork building step (seen at least 3 times in 10% of samples)
rpoc.usa <- filter_taxa(glom.rpoc.usa, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
# USA its
its.usa <- subset_samples(its.dat, Geog_loc == "USA")
glom.its.usa <- tax_glom(its.usa, "V8")
# filter low abundance taxa to simplify net.usawork building step (seen at least 3 times in 10% of samples)
its.usa<- filter_taxa(glom.its.usa, function(x) sum(x > 3) > (0.05*length(x)), TRUE)

#Niger rpoc
rpoc.nigeria <- subset_samples(rpoc.dat, Geog_loc == "Nigeria")
glom.rpoc.nigeria <- tax_glom(rpoc.nigeria, "V8")
# filter low abundance taxa to simplify net.usawork building step (seen at least 3 times in 10% of samples)
rpoc.nigeria <- filter_taxa(glom.rpoc.nigeria, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
# Niger ITS
its.nigeria <- subset_samples(its.dat, Geog_loc == "Nigeria")
glom.its.nigeria <- tax_glom(its.nigeria, "V8")
# filter low abundance taxa to simplify net.usawork building step (seen at least 3 times in 10% of samples)
its.nigeria<- filter_taxa(glom.its.nigeria, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
#net.usawork generation
pargs <- list(rep.num = 100, ncores=7, seed=123456)


usa.cross <- spiec.easi(list(rpoc.usa, its.usa), method='mb', nlambda=80,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(usa.cross)

set.seed(23123123)
nigeria.cross <- spiec.easi(list(rpoc.nigeria, its.nigeria), method='mb', nlambda=78,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(nigeria.cross)

usa.bact <- spiec.easi(list(rpoc.usa), method='mb', nlambda=70,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(usa.bact)

nigeria.bact <- spiec.easi(list(rpoc.nigeria), method='mb', nlambda=75,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(nigeria.bact)

save.image("cross_net.RData")





















































load("cross_net.RData")
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
friendly_pal("muted_nine", 9, type = "discrete")

import_taxa = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa","Stereum_hirsutum",
	"Streptococcus_mutans", "Treponema_denticola", 
 "Corynebacterium_durum","Leptotrichia_sp._oral_taxon_212", "Leptotrichia_sp._oral_taxon_215", "Corynebacterium_matruchotii", "Streptococcus_sanguinis")
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
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Leptotrichia_sp._oral_taxon_212"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "red"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Leptotrichia_sp._oral_taxon_215"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "red"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Corynebacterium_matruchotii"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "#AA4499"
neigh.nodes <- unlist(ego(net.usa, V(net.usa)[genus=="Streptococcus_sanguinis"], mode="all", order = 1))
inc.edges <- unlist(incident_edges(net.usa, neigh.nodes, mode = "all"))
ecol[inc.edges] <- "green"

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
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Leptotrichia_sp._oral_taxon_212"], mode="all")
ecol[inc.edges] <- "red"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Leptotrichia_sp._oral_taxon_215"], mode="all")
ecol[inc.edges] <- "red"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Corynebacterium_matruchotii"], mode="all")
ecol[inc.edges] <- "#AA4499"
inc.edges <- incident(net.usa,  V(net.usa)[genus=="Streptococcus_sanguinis"], mode="all")
ecol[inc.edges] <- "green"
#edge width
import_taxa2 = c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa","Stereum_hirsutum",
	"Streptococcus_mutans","Corynebacterium_durum", "Corynebacterium_matruchotii", "Leptotrichia_sp._oral_taxon_212", "Leptotrichia_sp._oral_taxon_215",
  "Streptococcus_sanguinis")
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
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Leptotrichia_sp._oral_taxon_212"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("red"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Leptotrichia_sp._oral_taxon_215"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("red"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Corynebacterium_matruchotii"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("#AA4499"), alpha =.6)
neigh.nodes <- neighbors(net.usa, V(net.usa)[genus=="Streptococcus_sanguinis"], mode="all")
vcol[neigh.nodes] <- adjustcolor(c("green"), alpha =.6)

vcol[V(net.usa)$genus=="Debaryomyces_prosopidis"] <-  adjustcolor(c("#332288"), alpha =.9)
vcol[V(net.usa)$genus=="Candida_tropicalis"] <-  adjustcolor(c("#117733"), alpha =.9)
vcol[V(net.usa)$genus=="Candida_albicans"] <-  adjustcolor(c("#CC6677"), alpha =.9)
vcol[V(net.usa)$genus=="Candida_unknown"] <-  adjustcolor(c("#88CCEE"), alpha =.9)
vcol[V(net.usa)$genus=="Rhodotorula_mucilaginosa"] <-  adjustcolor(c("#999933"), alpha =.9)
vcol[V(net.usa)$genus=="Stereum_hirsutum"] <-  adjustcolor(c("#882255"), alpha =.9)
vcol[V(net.usa)$genus=="Streptococcus_mutans"] <-  adjustcolor(c("#44AA99"), alpha =.9)
vcol[V(net.usa)$genus=="Treponema_denticola"] <-  adjustcolor(c("#DDCC77"), alpha =.9)
vcol[V(net.usa)$genus=="Corynebacterium_durum"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.usa)$genus=="Leptotrichia_sp._oral_taxon_212"] <-  adjustcolor(c("red"), alpha =.9)
vcol[V(net.usa)$genus=="Leptotrichia_sp._oral_taxon_215"] <-  adjustcolor(c("red"), alpha =.9)
vcol[V(net.usa)$genus=="Corynebacterium_matruchotii"] <-  adjustcolor(c("#AA4499"), alpha =.9)
vcol[V(net.usa)$genus=="Streptococcus_sanguinis"] <-  adjustcolor(c("green"), alpha =.9)

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

#graph
pdf("usa_positive_test.pdf", height = 100, width = 100)
plot(net.usa,vertex.color=vcol, edge.color=ecol, edge.width=ew, vertex.label=V(net.usa)$genus, vertex.size=vs, vertex.frame.color =vbod,
 edge.arrow.mode=0, vertex.frame.color = "white", vertex.label.cex = font, vertex.label.color = fontc, vertex.label.font=fontt)
dev.off()

llec <- cluster_louvain(net.usa, weights =NA)
#usa.louvain <- llec #GOOOD ONE
df <- as.data.frame(membership(llec))
df$names <- rownames(df)
df %>%  filter(names == "rASV2" | names == "ASV4" | names == "ASV13" | names == "rASV400" |names ==  "rASV864" |names ==  "rASV117"| names == "rASV103" | names =="rASV1327")
df %>%  filter(x ==6)


wc <- cluster_leading_eigen(net.usa, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
df %>%  filter(names == "ASV13" | names == "ASV13" | names == "rASV400" |names ==  "rASV864" |names ==  "rASV117"| names == "rASV105" | names =="rASV1327")
df %>%  filter(x ==4) %>%left_join(genus, by=c('names'='ASV'))




nllec <- cluster_louvain(net.nigeria, weights =NA)
df <- as.data.frame(membership(nllec))
df$names <- rownames(df)
df %>%  filter(names == "ASV4" |  names == "rASV339" |names ==  "rASV1905" |names ==  "rASV117" | names == "rASV104" | names == "rASV244")
df %>%  filter(x ==4)

nwc <- cluster_leading_eigen(net.nigeria, weights =NA)
df <- as.data.frame(membership(nwc))
df$names <- rownames(df)
df %>%  filter(names == "ASV4" |  names == "rASV339" |names ==  "rASV1905" |names ==  "rASV117" | names == "rASV104" | names == "rASV244")
df %>%  filter(x ==3) %>%left_join(genus, by=c('names'='ASV'))


save.image("test_goodusa.RData")
pdf("test2.pdf", width =100, height =100)
plot(wc, net.usa, vertex.label=V(net.usa)$genus, vertex.size=vs)
dev.off()

pdf("test.pdf", width =100, height =100)
plot(nwc, net.nigeria, vertex.label=V(net.nigeria)$genus, vertex.size=2)
dev.off()
pdf("test1.pdf", width =100, height =100)
plot(wc, net.usa, vertex.label=V(net.usa)$genus,vertex.size=vs)
dev.off()
#stats
#degree dist

dd.usa <- degree.distribution(net.usa.usa)
dd.nigeria <- degree.distribution(net.usa.nigeria)

pdf("deg_des.pdf")
plot(0:(length(dd.usa)-1), dd.usa, ylim=c(0,.35), type='b',
      ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.nigeria)-1), dd.nigeria, col="red" , type='b')
legend("topright", c("USA", "Nigeria"),
        col=c("black", "red"), pch=1, lty=1)
dev.off()


bact <- V(net.usa)[domain == "Bacteria"]$name
fungi <- V(net.usa)[domain == "Fungi"]$name

E(net.usa)[bact %--% fungi]
E(net.usa)[bact %--% bact]
E(net.usa)[fungi %--% fungi]


betweenness(
  net.usa,
  v = V(net.usa)$domain,
  directed = TRUE,
  weights = NULL,
  nobigint = TRUE,
  normalized = FALSE,
  cutoff = -1
)

etest <- rep("all", ecount(net.usa))
inc.edges <- incident(net.usa,  V(net.usa)[domain=="Bacteria"], mode="all")
etest[inc.edges] <- "bact"



