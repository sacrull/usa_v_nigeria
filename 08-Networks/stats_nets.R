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
setwd("/home/suzanne/usa_nigeria/networks/")
load("/home/suzanne/usa_nigeria/networks/cross_net.RData")
#stats
#degress dist
dd.usa <- degree.distribution(net.usa)
dd.nigeria <- degree.distribution(net.nigeria)
pdf("deg_des.pdf")
plot(0:(length(dd.usa)-1), dd.usa, ylim=c(0,.35), type='b',
      ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.nigeria)-1), dd.nigeria, col="red" , type='b')
legend("topright", c("USA", "Nigeria"),
        col=c("black", "red"), pch=1, lty=1)
dev.off()

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
#robustness
#betweeness
robust.usa <- robustness(net.usa, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
robust.usa.bact <- robustness(net.usa.bact, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
robust.nigeria <- robustness(net.nigeria, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
robust.nigeria.bact <- robustness(net.nigeria.bact, type = c("vertex"), measure = c("btwn.cent"), N = 1000)
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
  scale_y_continuous(name="Ratio of  remaining maximal component size to  initial maximal component size")+
  theme_minimal()
dev.off()
#degree
robust.usa <- robustness(net.usa, type = c("vertex"), measure = c("degree"), N = 1000)
robust.usa.bact <- robustness(net.usa.bact, type = c("vertex"), measure = c("degree"), N = 1000)
robust.nigeria <- robustness(net.nigeria, type = c("vertex"), measure = c("degree"), N = 1000)
robust.nigeria.bact <- robustness(net.nigeria.bact, type = c("vertex"), measure = c("degree"), N = 1000)
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
geom_label_repel(aes(label = Name),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
pdf("usa_nodes_btwn.pdf", width =20, height = 20)
ggplot(usa.comb, aes(x = nodes, y = between)) +
  geom_jitter(aes(color = phylum, shape = kingdom, size = 2), position = position_jitter(seed = 1))+
  scale_color_manual(values = mycolors) +
  geom_text(label=usa.comb$species, check_overlap = F, size = .2, position = position_jitter(seed = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  geom_xsidedensity(
    aes(
      y    = after_stat(density),
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
      x    = after_stat(density),
      yfill = factor(kingdom, levels =ords)
    ),
    alpha    = 0.5,
    size     = 1,
   # position = "stack"
  )+
  scale_yfill_manual(values = coul)+
  theme_bw()
dev.off()
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

pdf("nigeria_nodes_btwn.pdf", width =20, height = 20)
ggplot(nigeria.comb, aes(x = nodes, y = between)) +
  geom_jitter(aes(color = phylum, shape = kingdom, size = 2), position = position_jitter(seed = 1))+
  scale_color_manual(values = mycolors) +
  geom_text(label=nigeria.comb$species, check_overlap = F, size = .19, position = position_jitter(seed = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  geom_xsidedensity(
    aes(
      y    = after_stat(density),
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
      x    = after_stat(density),
      yfill = factor(kingdom, levels =ords)
    ),
    alpha    = 0.5,
    size     = 1,
    add = "mean"
    #position = "stack"
  )+
  scale_yfill_manual(values = coul)+
  theme_bw()
dev.off()
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
#bipart graph
bact <- V(net.usa)[domain == "Bacteria"]$name
fungi <- V(net.usa)[domain == "Fungi"]$name
fung.bact.edges <- E(net.usa)[bact %--% fungi]
fung.bact <- subgraph.edges(net.usa, eids = fung.bact.edges, delete.vertices = TRUE)
#make graph
V(fung.bact)$type <- ifelse(V(fung.bact)$domain == "Fungi", TRUE, FALSE)
is_bipartite(fung.bact)
V(fung.bact)$color <- ifelse(V(fung.bact)$type == TRUE, "orange", "cyan4")
pdf("bipart_usa.pdf")
plot(fung.bact, layout=layout_as_bipartite, arrow.mode=0, vertex.label=V(fung.bact)$genus, vertex.label.cex = 0.05, vertex.size=2, asp=0.5)
dev.off()
#highest node degree
degree.usa <- as.data.frame(degree(fung.bact, mode="all")) 
colnames(degree.usa)[1] <- "degree"
degree.usa$vector <- dplyr::pull(degree.usa, degree)
degree.usa[order(as.numeric(degree.usa[,2])), ]

#get only fungus and bact connections
bact <- V(net.nigeria)[domain == "Bacteria"]$name
fungi <- V(net.nigeria)[domain == "Fungi"]$name
fung.bact.edges <- E(net.nigeria)[bact %--% fungi]
fung.bact <- subgraph.edges(net.nigeria, eids = fung.bact.edges, delete.vertices = TRUE)
#make graph
V(fung.bact)$type <- ifelse(V(fung.bact)$domain == "Fungi", TRUE, FALSE)
is_bipartite(fung.bact)
V(fung.bact)$color <- ifelse(V(fung.bact)$type == TRUE, "orange", "cyan4")
pdf("bipart_nigeria.pdf", width=20, height=20)
plot(fung.bact, layout=layout_as_bipartite, arrow.mode=0, vertex.label=V(fung.bact)$genus, vertex.label.cex = 0.05, vertex.size=2, asp=0.5)
dev.off()
#highest node degree
degree.nigeria <- as.data.frame(degree(fung.bact, mode="all")) 
colnames(degree.nigeria)[1] <- "degree"
degree.nigeria$vector <- dplyr::pull(degree.nigeria, degree)
degree.nigeria[order(as.numeric(degree.nigeria[,2])), ]

#see if any same nodes
E(net.usa %s% net.nigeria)


#matrix connect network
bact <- V(net.usa)[domain == "Bacteria"]$name
fungi <- V(net.usa)[domain == "Fungi"]$name
fung.bact.edges <- E(net.usa)[bact %--% fungi]
fung.bact <- subgraph.edges(net.usa, eids = fung.bact.edges, delete.vertices = TRUE)
mat <- reshape2::melt(as.matrix(as_adjacency_matrix(fung.bact)))

mat <- filter(mat, value > 0)
#rearrange so all bact are in one column and fungi in other
mat1 <- mat[grepl("rASV", mat$Var1),] #get rows with rASV in Var1 sine it is mirrored on the other side
mat1 <- mat[grep("rASV", mat$Var1),]
pdf("cross_conects.usa.pdf")
ggplot(mat1, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      #coord_fixed()+
      theme_bw() +
      scale_fill_gradient(low="white", high="purple") +
      labs(y= "Fungi", x = "Bacteria")+
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        aspect.ratio = 1, # Force the plot into a square aspect ratio
        legend.position = "none") # Hide the legend (optional)
dev.off()

bact <- V(net.nigeria)[domain == "Bacteria"]$name
fungi <- V(net.nigeria)[domain == "Fungi"]$name
fung.bact.edges <- E(net.nigeria)[bact %--% fungi]
fung.bact <- subgraph.edges(net.nigeria, eids = fung.bact.edges, delete.vertices = TRUE)
mat <- reshape2::melt(as.matrix(as_adjacency_matrix(fung.bact)))

mat <- filter(mat, value > 0)
#rearrange so all bact are in one column and fungi in other
mat1 <- mat[grep("rASV", mat$Var1),]
pdf("cross_conects.nigeria.pdf")
ggplot(mat1, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      #coord_fixed()+
      theme_bw() +
      scale_fill_gradient(low="white", high="purple") +
      labs(y= "Fungi", x = "Bacteria")+
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        aspect.ratio = 1, # Force the plot into a square aspect ratio
        legend.position = "none") # Hide the legend (optional)
dev.off()
#connectance
connectance <- ecount(net.usa) / vcount(net.usa)^2 
connectance <- ecount(net.usa.bact) / vcount(net.usa.bact)^2 
connectance <- ecount(net.nigeria) / vcount(net.nigeria)^2 
connectance <- ecount(net.nigeria.bact) / vcount(net.nigeria.bact)^2 

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

#nestdness
adj.usa <- as.matrix(as_adjacency_matrix(net.usa))
adj.usa <- adj.usa[rowSums(adj.usa[,-1]) != 0,] #remove rows sum of 0
adj.usa <- adj.usa[, which(colSums(adj.usa) != 0)] #remove cols sum of 0

maxnodf(adj.usa, quality = 2)
