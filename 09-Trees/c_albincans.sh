grep -w "Candida_albicans" ../its/its_taxonomy.txt | awk '{print $1}' | sed 's/"//g' > candida_albicans.ids
seqtk subseq ../its/rep_set_its.fa candida_albicans.ids > candida_albicans.fa
mafft --auto candida_albicans.fa > candida_albicans.align.fa
raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100  -b 100 -x 02938 -n candida_albicans.tre -s candida_albicans.align.fa
iqtree -s candida_albicans.align.fa -m MFP -b 100 -nt 7


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


tre <- read.tree("candida_albicans.root.tre")

tre_dat <- merge_phyloseq(tre, rpoc.pd)
pdf("candida_albicans.pdf")
plot_tree(tre_dat, label.tips="taxa_names", color="Geog_loc", ladderize="left")
dev.off()

spec.dat <- subset_taxa(rpoc.pd, V8=="Candida_albicans")
tip_order <- rev(tre$tip.label) #tip order

data <- psmelt(spec.dat) # create dataframe from phyloseq object

pdf("./asv_dis.candida_albicans.pdf")
ggplot(data,aes(x=factor(OTU),y=Abundance,fill=factor(Geog_loc))) + 
  geom_bar(position="fill", stat="identity") + 
  scale_x_discrete(limits = tip_order,guide = guide_axis(angle = 90))+
  theme_minimal()
dev.off()

