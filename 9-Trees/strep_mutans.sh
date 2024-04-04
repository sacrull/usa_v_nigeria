grep -w "Streptococcus_mutans" ../rpoc/taxonomy_bac.txt | awk '{print $1}' | sed 's/"//g' > strep_mutans.ids
seqtk subseq ../rpoc/rep_set_rpoc.fa strep_mutans.ids > strep_mutans.fa
mafft --auto strep_mutans.fa > strep_mutans.align.fa
raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n strep_mutans.tre -s strep_mutans.align.fa

library(devtools)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ape)
library(microbiome)
library(ggpubr)
library(rstatix)

set.seed(12349)
setwd("~/usa_nigeria/tree")
load("~/usa_nigeria/phyloseq_obj/ps.RData")


tre <- read.tree("step_mutans.root.tre")

tre_dat <- merge_phyloseq(tre, rpoc.pd)
pdf("strep_mutans.pdf")
plot_tree(tre_dat, label.tips="taxa_names", color="Geog_loc", ladderize="left")
dev.off()

spec.dat <- subset_taxa(rpoc.pd, V8=="Streptococcus_mutans")
tip_order <- rev(tre$tip.label) #tip order

data <- psmelt(spec.dat) # create dataframe from phyloseq object

pdf("./asv_dis.lepto_212.pdf")
ggplot(data,aes(x=factor(OTU),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = tip_order,guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()