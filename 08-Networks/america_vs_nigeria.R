setwd("/home/suzanne/usa_nigeria/networks")
#load libraries
#devtools::install_github("stefpeschel/NetCoMi", 
#                         dependencies = c("Depends", "Imports", "LinkingTo"),
#                         repos = c("https://cloud.r-project.org/",
#                                   BiocManager::repositories()))
library(phyloseq)
library(NetCoMi)
library(tidyverse)
library(reshape2)
library(limma)
#library(ciriclize)
#load data
load("../phyloseq_obj/ps.RData")
set.seed=10010
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
pargs <- list(rep.num = 100, ncores=7, seed=10010)


usa.cross <- spiec.easi(list(rpoc.usa, its.usa), method='mb', nlambda=40,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))

nigeria.cross <- spiec.easi(list(rpoc.nigeria, its.nigeria), method='mb', nlambda=40,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))

save.image("cross_nets.RData")

library(SpiecEasi)
library(igraph)
library(viridis)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(Matrix)
load("cross_nets.RData")

#stats
#statistics
getStability(usa.cross)
getStability(nigeria.cross)
#usa.cross
sebeta.usa.cross <- symBeta(getOptBeta(usa.cross), mode='maxabs')
taxnames <- c(taxa_names(rpoc.usa), taxa_names(its.usa))
colnames(sebeta.usa.cross) <- rownames(sebeta.usa.cross) <- taxnames
weighted.adj.usa.cross <- sebeta.usa.cross*getRefit(usa.cross)
#run weighted object
grph.usa <- adj2igraph(weighted.adj.usa.cross)
#nigeria.cross
sebeta.nigeria.cross <- symBeta(getOptBeta(nigeria.cross), mode='maxabs')
taxnames <- c(taxa_names(rpoc.nigeria), taxa_names(its.nigeria))
colnames(sebeta.nigeria.cross) <- rownames(sebeta.nigeria.cross) <- taxnames
weighted.adj.nigeria.cross <- sebeta.nigeria.cross*getRefit(nigeria.cross)
#run weighted object
grph.nigeria <- adj2igraph(weighted.adj.nigeria.cross)

library(circlize)
library(tidyverse)
library(reshape2)
#getting just genus name
taxa = as(tax_table(ps.data), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
orderdf <- orderdf %>% 
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

#replace asv with spiecieas
spiece_dataframe1.hcf <- left_join(adjacencyData.usa, orderdf, by=c('from'='ASV'))
colnames(spiece_dataframe1.hcf)[4] <- "origin"
spiece_dataframe2.hcf <- left_join(spiece_dataframe1.hcf, orderdf, by=c('to'='ASV'))
colnames(spiece_dataframe2.hcf)[5] <- "destination"
spiece_dataframe3.hcf <- select(spiece_dataframe2.hcf, weight, origin, destination)
colnames(spiece_dataframe3.hcf)[1] <- "value"
spiece_dataframe3.2.hcf <- spiece_dataframe3.hcf[, c(2,3,1)]
spiece_dataframe3.3.hcf <- subset(spiece_dataframe3.2.hcf, value != 0) #remove no connections

write.table(spiece_dataframe3.3.hcf, "usa_cross.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


pdf("corss_usa.pdf")
chordDiagram(spiece_dataframe3.2.hcf, order = ordered_list2, transparency = 0.5, col = ifelse(spiece_dataframe3.2.hcf$origin == "Candida_albicans" | spiece_dataframe3.2.hcf$origin == "Leptotrichia_sp._oral_taxon_212" | spiece_dataframe3.2.hcf$origin =="Leptotrichia_sp._oral_taxon_215" | spiece_dataframe3.2.hcf$destination == "Candida_albicans" | spiece_dataframe3.2.hcf$destination == "Leptotrichia_sp._oral_taxon_212" | spiece_dataframe3.2.hcf$destination =="Leptotrichia_sp._oral_taxon_215" | spiece_dataframe3.2.hcf$origin== "Streptococcus_mutans" | spiece_dataframe3.2.hcf$destination== "Streptococcus_mutans", ifelse(spiece_dataframe3.2.hcf$origin== "Streptococcus_mutans" | spiece_dataframe3.2.hcf$destination== "Streptococcus_mutans", "red", "green"), "gray"),
  annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) < 20) {
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
  niceFacing = TRUE, adj = c(0, 0.5), cex=.3)
  } else {
  circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
  niceFacing = TRUE, adj = c(0.5, 0), cex=.3)
  }
}, bg.border = NA)
dev.off()

#making data frame from igraph
#nigeria
spiec_df_usa <- as.matrix(weighted.adj.nigeria.cross)
adjacencyData.usa <- as.data.frame(melt(spiec_df_usa))
colnames(adjacencyData.usa)[1] <- "to"
colnames(adjacencyData.usa)[2] <- "from"
colnames(adjacencyData.usa)[3] <- "weight"

#replace asv with spiecieas
spiece_dataframe1.hcf <- left_join(adjacencyData.usa, orderdf, by=c('from'='ASV'))
colnames(spiece_dataframe1.hcf)[4] <- "origin"
spiece_dataframe2.hcf <- left_join(spiece_dataframe1.hcf, orderdf, by=c('to'='ASV'))
colnames(spiece_dataframe2.hcf)[5] <- "destination"
spiece_dataframe3.hcf <- select(spiece_dataframe2.hcf, weight, origin, destination)
colnames(spiece_dataframe3.hcf)[1] <- "value"
spiece_dataframe3.2.hcf <- spiece_dataframe3.hcf[, c(2,3,1)]
spiece_dataframe3.3.hcf <- subset(spiece_dataframe3.2.hcf, value != 0) #remove no connections

write.table(spiece_dataframe3.3.hcf, "nigeria_cross.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

pdf("corss_nigeria.pdf")
chordDiagram(spiece_dataframe3.3.hcf, order = ordered_list2, transparency = 0.5, col = ifelse(spiece_dataframe3.2.hcf$origin == "Candida_albicans" | spiece_dataframe3.2.hcf$origin == "Leptotrichia_sp._oral_taxon_212" | spiece_dataframe3.2.hcf$origin =="Leptotrichia_sp._oral_taxon_215" | spiece_dataframe3.2.hcf$destination == "Candida_albicans" | spiece_dataframe3.2.hcf$destination == "Leptotrichia_sp._oral_taxon_212" | spiece_dataframe3.2.hcf$destination =="Leptotrichia_sp._oral_taxon_215" | spiece_dataframe3.2.hcf$origin== "Streptococcus_mutans" | spiece_dataframe3.2.hcf$destination== "Streptococcus_mutans", ifelse(spiece_dataframe3.2.hcf$origin== "Streptococcus_mutans" | spiece_dataframe3.2.hcf$destination== "Streptococcus_mutans", "red", "green"), "gray"),
  annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) < 20) {
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
  niceFacing = TRUE, adj = c(0, 0.5), cex=.3)
  } else {
  circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
  niceFacing = TRUE, adj = c(0.5, 0), cex=.3)
  }
}, bg.border = NA)
dev.off()


#anuran
sed '1d' nigeria_cross.txt > ./geo_nets/nigeria_cross.txt
sed '1d' usa_cross.txt > ./geo_nets/usa_cross.txt
conda activate anuran
anuran -i geo_nets -o anuran_out -size 0.6 1 -perm 5 -nperm 10 -draw
# also make network with all samples -- size minimum 50% of the full gropus
anuran -i all_networks -o anuran_ALL_out -draw -perm 5 -nperm 10 -size 0.5 -c









as.data.frame(tax_table(rpoc.usa))$V8