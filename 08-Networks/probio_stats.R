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
library(smplot2)
library(foreach)
setwd("/home/suzanne/usa_nigeria/networks/probio")
load("/home/suzanne/usa_nigeria/networks/probio/cross_net.RData")
#stats
#degress dist
dd.usa <- degree.distribution(net.usa)
dd.nigeria <- degree.distribution(net.nigeria)
dd.usa.bact <- degree.distribution(net.usa.bact)
dd.nigeria.bact <- degree.distribution(net.nigeria.bact)
pdf("deg_des.pdf")
plot(0:(length(dd.usa)-1), dd.usa, col = "#6EBDC2", ylim=c(0,.35), type='b',
      ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.nigeria)-1), dd.nigeria, col="#DC756D" , type='b')
points(0:(length(dd.usa.bact)-1), dd.usa.bact, col="lightblue" , type='b')
points(0:(length(dd.nigeria.bact)-1), dd.nigeria.bact, col="darkorange" , type='b')
legend("topright", c("USA_Cross", "Nigeria_Cross", "USA_Bact", "Nigeria_Bact"),
        col=c("#6EBDC2", "#DC756D", "lightblue", "darkorange"), pch=1, lty=1)
dev.off()

#clusters
#wc <- cluster_leading_eigen(net.usa, weights =NA)
test.net <- delete.edges(net.usa, E(net.usa)[E(net.usa)$weight >0.3])
wc <- cluster_fast_greedy(test.net, weights =NA)
df <- as.data.frame(membership(wc))
df$names <- rownames(df)
#df %>%  filter(names == "rASV2" | names == "ASV4" | names == "rASV400" |names ==  "rASV864" |names ==  "rASV117"| names == "rASV105" | names =="rASV1327")
#df %>%  filter(x ==5) %>%left_join(genus, by=c('names'='ASV'))
#df %>%left_join(genus, by=c('names'='ASV')) %>% filter(V8 == "Streptococcus_parasanguinis"| V8 == "Limosilactobacillus_fermentum" | V8 == "Streptococcus_mutans" | 
#	V8 == "Scardovia_wiggsiae"| V8 == "Rhodotorula_mucilaginosa" | V8 == "Debaryomyces_prosopidis" | V8 =="Leptotrichia_sp._oral_taxon_212")

df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | 
  V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum")

df %>%  filter(x ==1) %>%left_join(genus, by=c('names'='ASV'))

nwc <- cluster_leading_eigen(net.nigeria, weights =NA)
test.net <- delete.edges(net.nigeria, E(net.nigeria)[E(net.nigeria)$weight >0.3])
nwc <- cluster_fast_greedy(test.net, weights =NA)
df <- as.data.frame(membership(nwc))
df$names <- rownames(df)
#df %>%  filter(names == "rASV3" |  names == "rASV339" |names ==  "rASV1905" |names ==  "rASV117" | names == "rASV104" | names == "rASV244")
#df %>%  filter(x ==6) %>%left_join(genus, by=c('names'='ASV'))
df %>%left_join(genus, by=c('names'='ASV')) %>% filter( V8 == "Debaryomyces_prosopidis" | V8 == "Streptococcus_sanguinis"| V8 == "Streptococcus_oralis" | 
  V8 == "Streptococcus_parasanguinis" | V8 =="Leptotrichia_sp._oral_taxon_212" | V8 =="Leptotrichia_sp._oral_taxon_215" | V8 == "Corynebacterium_durum")

df %>%  filter(x ==6) %>%left_join(genus, by=c('names'='ASV'))

#connections between domains
#usa
bact <- V(net.usa)[domain == "Bacteria"]$name #label bacteria
fungi <- V(net.usa)[domain == "Fungi"]$name #label fungus
bact.bact <- gsize(subgraph.edges(net.usa, eids = E(net.usa)[bact %--% bact], delete.vertices = TRUE)) #count edges
fung.fung <- gsize(subgraph.edges(net.usa, eids = E(net.usa)[fungi %--% fungi], delete.vertices = TRUE))
fung.bact <- gsize(subgraph.edges(net.usa, eids = E(net.usa)[bact %--% fungi], delete.vertices = TRUE))
edges.usa <- data.frame(bact.bact, fung.fung, fung.bact) #make dataframes
#nigeria
bact <- V(net.nigeria)[domain == "Bacteria"]$name
fungi <- V(net.nigeria)[domain == "Fungi"]$name
bact.bact <- gsize(subgraph.edges(net.nigeria, eids = E(net.nigeria)[bact %--% bact], delete.vertices = TRUE))
fung.fung <- gsize(subgraph.edges(net.nigeria, eids = E(net.nigeria)[fungi %--% fungi], delete.vertices = TRUE))
fung.bact<- gsize(subgraph.edges(net.nigeria, eids = E(net.nigeria)[bact %--% fungi], delete.vertices = TRUE))
edges.nigeria <- data.frame(bact.bact, fung.fung, fung.bact) #make dataframes
#combine
edges.usa <- mutate(edges.usa, location = "USA") #adding usa column
edges.nigeria <- mutate(edges.nigeria, location = "Nigeria") #adding rpoc column
edges_counts <- melt(rbind(edges.usa, edges.nigeria)) #combining and melting
pdf("./net_edges_barchart.pdf")
ggplot(edges_counts, aes(fill=variable, y=value, x=location)) + geom_bar(position="fill", stat="identity") + theme_minimal()
dev.off()

#nodes vs betweeness
#america
#nodes
nodes.usa <- igraph::degree(net.usa, v = V(net.usa), mode = c("all"), loops = TRUE, normalized = FALSE)
nodes.usa <- as.data.frame(nodes.usa) 
nodes.usa <- nodes.usa %>% 
  rownames_to_column(var = "ASV")
colnames(nodes.usa)[2] <- "nodes"
#between
btwn.usa <- betweenness(net.usa,v = V(net.usa),directed = FALSE, weights = NA, normalized = FALSE, cutoff = -1)
#eig.usa <- eigen_centrality(net.usa, directed = FALSE, scale = TRUE,
#  weights = NULL, options = arpack_defaults)
#btwn.usa <- as.data.frame(eig.usa$vector)
btwn.usa <- as.data.frame(btwn.usa) 
btwn.usa <- btwn.usa %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.usa)[2] <- "between"
#combine
usa.comb <- left_join(btwn.usa, nodes.usa, by=c('ASV'='ASV'))
#rename to domain (domain was created in a previous session check usa_net_vis.R)
usa.comb <- left_join(usa.comb, domain, by=c('ASV'='ASV'))
colnames(usa.comb)[4] <- "kingdom"
#getting phyla names
taxa = as(tax_table(ps.data), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V3)
phyla <- orderdf %>% 
  rownames_to_column(var = "ASV")
#rename to phyla
usa.comb <- left_join(usa.comb, phyla, by=c('ASV'='ASV'))
colnames(usa.comb)[5] <- "phylum"
usa.comb <- filter(usa.comb, between > 0)
#rename to species
usa.comb <- left_join(usa.comb, genus, by=c('ASV'='ASV'))
colnames(usa.comb)[6] <- "species"
#colors
nb.cols <- 13
mycolors <-  colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
coul <-c("orange", "cyan4")
ords <- c("Fungi", "Bacteria")
usa_lavs <- c("Debaryomyces_prosopidis", "Candida_tropicalis", "Candida_albicans", "Candida_unknown","Rhodotorula_mucilaginosa",
	"Stereum_hirsutum","Streptococcus_mutans", "Treponema_denticola", "Corynebacterium_durum",
	"Corynebacterium_segmentosum", "Corynebacterium_mustelae", "Corynebacterium_matruchotii")
#count by kingdom
usa.comb %>% group_by(kingdom) %>% 
  summarise(total_count=n(),
            .groups = 'drop')
#nigeria
nodes.nigeria <- igraph::degree(net.nigeria, v = V(net.nigeria), mode = c("all"), loops = TRUE, normalized = FALSE)
nodes.nigeria <- as.data.frame(nodes.nigeria) 
nodes.nigeria <- nodes.nigeria %>% 
  rownames_to_column(var = "ASV")
colnames(nodes.nigeria)[2] <- "nodes"
#between
btwn.nigeria <- betweenness(net.nigeria,v = V(net.nigeria),directed = FALSE, weights = NA, normalized = FALSE, cutoff = -1)
btwn.nigeria <- as.data.frame(btwn.nigeria) 
btwn.nigeria <- btwn.nigeria %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.nigeria)[2] <- "between"
#combine
nigeria.comb <- left_join(btwn.nigeria, nodes.nigeria, by=c('ASV'='ASV'))
#rename to domain (domain was created in a previous session check nigeria_net_vis.R)
nigeria.comb <- left_join(nigeria.comb, domain, by=c('ASV'='ASV'))
colnames(nigeria.comb)[4] <- "kingdom"
#rename to phyla
nigeria.comb <- left_join(nigeria.comb, phyla, by=c('ASV'='ASV'))
colnames(nigeria.comb)[5] <- "phylum"
nigeria.comb <- filter(nigeria.comb, between > 0)
#rename to species
nigeria.comb <- left_join(nigeria.comb, genus, by=c('ASV'='ASV'))
colnames(nigeria.comb)[6] <- "species"

nigeria.comb %>% group_by(kingdom) %>% 
  summarise(total_count=n(),
            .groups = 'drop')
#median diff for nodes and betweeness
usa.comb <- mutate(usa.comb, location = "usa")
nigeria.comb <- mutate(nigeria.comb, location = "nigeria")
btwn.nets <- rbind(usa.comb, nigeria.comb)
btwn.nets <- btwn.nets %>%  na.omit()
pdf("between_diff.pdf")
ggplot(btwn.nets, aes(x=location,y=between))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_violin(aes(fill=location),alpha    = 0.5) +
  geom_jitter(aes(color=location), shape=16, position=position_jitter(0.2), size=2.5)+
  geom_boxplot(width=0.1)+
  labs(x ="Geographic Location", y = "Betweeness Centrality")+
  facet_wrap(~kingdom)+
  theme_classic()
dev.off()
sm_auc(x, subj_b_day1$Cbratio)

pdf("node_diff.pdf")
ggplot(btwn.nets, aes(x=location,y=nodes))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_violin(aes(fill=location),alpha    = 0.5) +
  geom_jitter(aes(color=location), shape=16, position=position_jitter(0.2), size=2.5)+
  geom_boxplot(width=0.1)+
  labs(x ="Geographic Location", y = "Betweeness Centrality")+
  facet_wrap(~kingdom)+
  theme_classic()
dev.off()
#histograms by total density
pdf("usa_nodes_btwn.pdf", width =20, height = 20)
ggplot(usa.comb, aes(x = nodes, y = between)) +
  geom_jitter(aes(color = phylum, shape = kingdom, size = 2), position = position_jitter(seed = 1))+
  scale_color_manual(values = mycolors) +
  geom_text(label=usa.comb$species, check_overlap = F, size = .2, position = position_jitter(seed = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  geom_xsidedensity(
    aes(
      y    = ..count../sum(..count..),
      xfill = factor(kingdom, levels =ords),
    ),
    alpha    = 0.5,
    size     = 1,
    #position = "stack"
  ) +
  scale_xsidey_continuous(minor_breaks = NULL)+
  scale_xfill_manual(values = coul)+
  geom_ysidedensity(
    aes(
      x    = ..count../sum(..count..),
      yfill = factor(kingdom, levels =ords)
    ),
    alpha    = 0.5,
    size     = 1,
   # position = "stack"
  )+
  scale_yfill_manual(values = coul)+
  theme_bw()
dev.off()
pdf("nigeria_nodes_btwn.pdf", width =20, height = 20)
ggplot(nigeria.comb, aes(x = nodes, y = between)) +
  geom_jitter(aes(color = phylum, shape = kingdom, size = 2), position = position_jitter(seed = 1))+
  scale_color_manual(values = mycolors) +
  geom_text(label=nigeria.comb$species, check_overlap = F, size = .19, position = position_jitter(seed = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  geom_xsidedensity(
    aes(
      y    = ..count../sum(..count..),
      xfill = factor(kingdom, levels =ords),
    ),
    alpha    = 0.5,
    size     = 1,
    #position = "stack"
  ) +
  scale_xsidey_continuous(minor_breaks = NULL)+
  scale_xfill_manual(values = coul)+
  geom_ysidedensity(
    aes(
      x    = ..count../sum(..count..),
      yfill = factor(kingdom, levels =ords)
    ),
    alpha    = 0.5,
    size     = 1,
    #position = "stack"
  )+
  scale_yfill_manual(values = coul)+
  theme_bw()
dev.off()

#modularity
wtc <- cluster_fast_greedy(net.usa, weights = NA)
modularity(net.usa, membership(wtc))
wtc <- cluster_fast_greedy(net.usa.bact, weights = NA)
modularity(net.usa.bact, membership(wtc))
wtc <- cluster_fast_greedy(net.nigeria, weights = NA)
modularity(net.nigeria, membership(wtc))
wtc <- cluster_fast_greedy(net.nigeria.bact, weights = NA)
modularity(net.nigeria.bact, membership(wtc))

#transitivity
transitivity(net.usa, type = c("global"), weights = NA)
transitivity(net.usa.bact, type = c("global"), weights = NA)
transitivity(net.nigeria, type = c("global"), weights = NA)
transitivity(net.nigeria.bact, type = c("global"), weights = NA)

#edge denisty
edge_density(net.usa, loops=TRUE) 
edge_density(net.usa.bact, loops=TRUE) 
edge_density(net.nigeria, loops=TRUE) 
edge_density(net.nigeria.bact, loops=TRUE) 


#betweeness centrality
btwn.usa <- betweenness(net.usa,v = V(net.usa),directed = FALSE, weights = NA, normalized = TRUE, cutoff = -1)
btwn.usa <- as.data.frame(btwn.usa) 
btwn.usa <- btwn.usa %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.usa)[2] <- "cross_btwn"
btwn.usa.bact <- betweenness(net.usa.bact,v = V(net.usa.bact),directed = FALSE, weights = NA, normalized = TRUE, cutoff = -1)
btwn.usa.bact <- as.data.frame(btwn.usa.bact) 
btwn.usa.bact <- btwn.usa.bact %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.usa.bact)[2] <- "bact_btwn"
#join the different networks
btwn.usa <- left_join(btwn.usa, btwn.usa.bact, by=c('ASV'='ASV'))
btwn.usa <- mutate(btwn.usa, location = "usa")
#nigeria
btwn.nigeria <- betweenness(net.nigeria,v = V(net.nigeria),directed = FALSE, weights = NA, normalized = TRUE, cutoff = -1)
btwn.nigeria <- as.data.frame(btwn.nigeria) 
btwn.nigeria <- btwn.nigeria %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.nigeria)[2] <- "cross_btwn"
btwn.nigeria.bact <- betweenness(net.nigeria.bact,v = V(net.nigeria.bact),directed = FALSE, weights = NA, normalized = TRUE, cutoff = -1)
btwn.nigeria.bact <- as.data.frame(btwn.nigeria.bact) 
btwn.nigeria.bact <- btwn.nigeria.bact %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.nigeria.bact)[2] <- "bact_btwn"
#join the different networks
btwn.nigeria <- left_join(btwn.nigeria, btwn.nigeria.bact, by=c('ASV'='ASV'))
btwn.nigeria <- mutate(btwn.nigeria, location = "nigeria")
#combine dataframes 
btwn.nets <- rbind(btwn.usa, btwn.nigeria)
btwn.nets <- btwn.nets %>%  na.omit()
#plot
pdf("betweeness.pdf")
ggscatter(btwn.nets, x = "bact_btwn", y = "cross_btwn",
   color = "location",
   add = "reg.line", conf.int = TRUE)+
 stat_cor(aes(color = location))
dev.off()

#eigen centrality
#usa
eig.usa <- eigen_centrality(net.usa, directed = FALSE, scale = TRUE,
  weights = NA, options = arpack_defaults)
btwn.usa <- as.data.frame(eig.usa$vector)
btwn.usa <- btwn.usa %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.usa)[2] <- "cross_eign"
eig.usa <- eigen_centrality(net.usa.bact, directed = FALSE, scale = TRUE,
  weights = NA, options = arpack_defaults)
btwn.usa.bact <- as.data.frame(eig.usa$vector)
btwn.usa.bact <- btwn.usa.bact %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.usa.bact)[2] <- "bact_eign"
#join the different networks
btwn.usa <- left_join(btwn.usa, btwn.usa.bact, by=c('ASV'='ASV'))
btwn.usa <- mutate(btwn.usa, location = "usa")
#nigeria
eig.nigeria <- eigen_centrality(net.nigeria, directed = FALSE, scale = TRUE,
  weights = NA, options = arpack_defaults)
btwn.nigeria <- as.data.frame(eig.nigeria$vector)
btwn.nigeria <- btwn.nigeria %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.nigeria)[2] <- "cross_eign"
eig.nigeria <- eigen_centrality(net.nigeria.bact, directed = FALSE, scale = TRUE,
  weights = NA, options = arpack_defaults)
btwn.nigeria.bact <- as.data.frame(eig.nigeria$vector)
btwn.nigeria.bact <- btwn.nigeria.bact %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.nigeria.bact)[2] <- "bact_eign"
#join the different networks
btwn.nigeria <- left_join(btwn.nigeria, btwn.nigeria.bact, by=c('ASV'='ASV'))
btwn.nigeria <- mutate(btwn.nigeria, location = "nigeria")
#combine data
btwn.nets <- rbind(btwn.usa, btwn.nigeria)
btwn.nets <- btwn.nets %>%  na.omit()
#plot
pdf("eigen.pdf")
ggscatter(btwn.nets, x = "bact_eign", y = "cross_eign",
   color = "location",
   add = "reg.line", conf.int = TRUE)+
 stat_cor(aes(color = location))
dev.off()

#clossness
#usa
btwn.usa <- closeness(net.usa, vids = V(net.usa), mode = c("all"), weights = NA, normalized = FALSE, cutoff = -1)
btwn.usa <- as.data.frame(btwn.usa) 
btwn.usa <- btwn.usa %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.usa)[2] <- "cross_close"
btwn.usa.bact <- closeness(net.usa.bact, vids = V(net.usa.bact), mode = c("all"), weights = NA, normalized = FALSE, cutoff = -1)
btwn.usa.bact <- as.data.frame(btwn.usa.bact) 
btwn.usa.bact <- btwn.usa.bact %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.usa.bact)[2] <- "bact_close"
#join the different networks
btwn.usa <- left_join(btwn.usa, btwn.usa.bact, by=c('ASV'='ASV'))
btwn.usa <- mutate(btwn.usa, location = "usa")
#nigeria
btwn.nigeria <- closeness(net.nigeria, vids = V(net.nigeria), mode = c("all"), weights = NA, normalized = FALSE, cutoff = -1)
btwn.nigeria <- as.data.frame(btwn.nigeria) 
btwn.nigeria <- btwn.nigeria %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.nigeria)[2] <- "cross_close"
btwn.nigeria.bact <- closeness(net.nigeria.bact, vids = V(net.nigeria.bact), mode = c("all"), weights = NA, normalized = FALSE, cutoff = -1)
btwn.nigeria.bact <- as.data.frame(btwn.nigeria.bact) 
btwn.nigeria.bact <- btwn.nigeria.bact %>% 
  rownames_to_column(var = "ASV")
colnames(btwn.nigeria.bact)[2] <- "bact_close"
#join the different networks
btwn.nigeria <- left_join(btwn.nigeria, btwn.nigeria.bact, by=c('ASV'='ASV'))
btwn.nigeria <- mutate(btwn.nigeria, location = "nigeria")
#combine dataframes 
btwn.nets <- rbind(btwn.usa, btwn.nigeria)
btwn.nets <- btwn.nets %>%  na.omit()
#plot
pdf("closseness.pdf")
ggscatter(btwn.nets, x = "bact_close", y = "cross_close",
   color = "location",
   add = "reg.line", conf.int = TRUE)+
 stat_cor(aes(color = location))
dev.off()

#robustness
#betweeness
robust.usa <- robustness(net.usa, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
sm_auc(robust.usa$removed.pct, robust.usa$comp.pct) #0.381739
robust.usa.bact <- robustness(net.usa.bact, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
sm_auc(robust.usa.bact$removed.pct, robust.usa.bact$comp.pct) #0.2682372
robust.nigeria <- robustness(net.nigeria, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
sm_auc(robust.nigeria$removed.pct, robust.nigeria$comp.pct) #0.3762612
robust.nigeria.bact <- robustness(net.nigeria.bact, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
sm_auc(robust.nigeria.bact$removed.pct, robust.nigeria.bact$comp.pct) #0.3642082
legend_colors <- c("Nigeria" = "#DC756D", "USA"= "#6EBDC2")
legend_line <- c("Cross_Kingdom"= "solid","Bacteria_only" = "dashed")
pdf("robust_btwn.pdf")
ggplot() +
  geom_line(data = robust.usa, aes(x=removed.pct, y=comp.pct, colour = "USA",linetype= "Cross_Kingdom", col=interaction(measure)), size =1.25) + # must include argument label "data"
  geom_line(data = robust.usa.bact, aes(x=removed.pct, y=comp.pct, colour = "USA", linetype= "Only_Bacteria", col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Cross_Kingdom",col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria.bact, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Only_Bacteria",col=interaction(measure)), size =1.25)+
  labs(colour = "Legend") + 
  scale_color_manual(values = legend_colors) +
  scale_x_continuous(name="Precentage of nodes removed") +
  scale_y_continuous(name="Ratio of remaining maximal component size to  initial maximal component size")+
  theme_minimal()
dev.off()
#degree
robust.usa <- robustness(net.usa, type = c("vertex"), measure = c("degree"), N = 1000)
sm_auc(robust.usa$removed.pct, robust.usa$comp.pct) #0.3737553
robust.usa.bact <- robustness(net.usa.bact, type = c("vertex"), measure = c("degree"), N = 1000)
sm_auc(robust.usa.bact$removed.pct, robust.usa.bact$comp.pct) #0.2562688
robust.nigeria <- robustness(net.nigeria, type = c("vertex"), measure = c("degree"), N = 1000)
sm_auc(robust.nigeria$removed.pct, robust.nigeria$comp.pct) #0.3571347
robust.nigeria.bact <- robustness(net.nigeria.bact, type = c("vertex"), measure = c("degree"), N = 1000)
sm_auc(robust.nigeria.bact$removed.pct, robust.nigeria.bact$comp.pct) #0.3379897
pdf("robust_degree.pdf")
ggplot() +
  geom_line(data = robust.usa, aes(x=removed.pct, y=comp.pct, colour = "USA",linetype= "Cross_Kingdom", col=interaction(measure)), size =1.25) + # must include argument label "data"
  geom_line(data = robust.usa.bact, aes(x=removed.pct, y=comp.pct, colour = "USA", linetype= "Only_Bacteria", col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Cross_Kingdom",col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria.bact, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Only_Bacteria",col=interaction(measure)), size =1.25)+
  labs(colour = "Legend") + 
  scale_color_manual(values = legend_colors) +
  scale_x_continuous(name="Precentage of nodes removed") +
  scale_y_continuous(name="Ratio of  remaining maximal component size to  initial maximal component size")+
  theme_minimal()
dev.off()
#random used updated function
robust.usa <- robustness(net.usa, type = c("vertex"), measure = c("random"), N = 1000)
robust.usa.bact <- robustness(net.usa.bact, type = c("vertex"), measure = c("random"), N = 1000)
robust.nigeria <- robustness(net.nigeria, type = c("vertex"), measure = c("random"), N = 1000)
robust.nigeria.bact <- robustness(net.nigeria.bact, type = c("vertex"), measure = c("random"), N = 1000)
pdf("robust_random.pdf")
ggplot() +
  geom_line(data = robust.usa, aes(x=removed.pct, y=comp.pct, colour = "USA",linetype= "Cross_Kingdom", col=interaction(measure)), size =1.25) + # must include argument label "data"
  geom_line(data = robust.usa.bact, aes(x=removed.pct, y=comp.pct, colour = "USA", linetype= "Only_Bacteria", col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Cross_Kingdom",col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria.bact, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Only_Bacteria",col=interaction(measure)), size =1.25)+
  labs(colour = "Legend") + 
  scale_color_manual(values = legend_colors) +
  scale_x_continuous(name="Precentage of nodes removed") +
  scale_y_continuous(name="Ratio of  remaining maximal component size to  initial maximal component size")+
  theme_minimal()
dev.off()
#eigen use custom function
robust.usa <- robustness_eigen(net.usa, type = c("vertex"), measure = c("eigen"), N = 1000)
robust.usa.bact <- robustness_eigen(net.usa.bact, type = c("vertex"), measure = c("eigen"), N = 1000)
robust.nigeria <- robustness_eigen(net.nigeria, type = c("vertex"), measure = c("eigen"), N = 1000)
robust.nigeria.bact <- robustness_eigen(net.nigeria.bact, type = c("vertex"), measure = c("eigen"), N = 1000)
pdf("robust_eigen.pdf")
ggplot() +
  geom_line(data = robust.usa, aes(x=removed.pct, y=comp.pct, colour = "USA",linetype= "Cross_Kingdom", col=interaction(measure)), size =1.25) + # must include argument label "data"
  geom_line(data = robust.usa.bact, aes(x=removed.pct, y=comp.pct, colour = "USA", linetype= "Only_Bacteria", col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Cross_Kingdom",col=interaction(measure)), size =1.25)+
  geom_line(data = robust.nigeria.bact, aes(x=removed.pct, y=comp.pct, colour = "Nigeria",linetype= "Only_Bacteria",col=interaction(measure)), size =1.25)+
  labs(colour = "Legend") + 
  scale_color_manual(values = legend_colors) +
  scale_x_continuous(name="Precentage of nodes removed") +
  scale_y_continuous(name="Ratio of  remaining maximal component size to  initial maximal component size")+
  theme_minimal()
dev.off()