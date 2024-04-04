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

#Niger rpoc
rpoc.nigeria <- subset_samples(rpoc.dat, Geog_loc == "Nigeria")
glom.rpoc.nigeria <- tax_glom(rpoc.nigeria, "V8")
# filter low abundance taxa to simplify network building step (seen at least 3 times in 10% of samples)
rpoc.nigeria <- filter_taxa(glom.rpoc.nigeria, function(x) sum(x > 3) > (0.05*length(x)), TRUE)

#network generation
pargs <- list(rep.num = 100, ncores=7, seed=123456)


usa.bact <- spiec.easi(list(rpoc.usa), method='mb', nlambda=70,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(usa.bact)

nigeria.bact <- spiec.easi(list(rpoc.nigeria), method='mb', nlambda=70,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=7))
getStability(nigeria.bact)

save.image("bact_nets.RData")