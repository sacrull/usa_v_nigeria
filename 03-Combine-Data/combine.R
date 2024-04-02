#library(devtools)
#BiocManager::install("phyloseq")
#BiocManager::install("Biostrings")

library(phyloseq)
library(Biostrings)
setwd("~/usa_nigeria/phyloseq_obj")
#its
#frequency
seqtab_its <- read.table("~/usa_nigeria/its/sequence_table_its.merged.txt", header=T, row.names=1)
asv_its <-otu_table(seqtab_its, taxa_are_rows=F)
#taxonomy
tax_tab_its <- read.table("~/usa_nigeria/its/its_taxonomy.txt", header=F, row.names=1, sep="\t")
tax_its <- tax_table(as.matrix(tax_tab_its))
#ref sequences
refseq_its <- Biostrings::readDNAStringSet("~/usa_nigeria/its/rep_set_its.fa")
#metadata
metadata <- read.table("~/usa_nigeria/Nigeria_USA_meta.txt", sep="\t", header=T, row.names=1)
map <- sample_data(metadata)
#merge into one phyloseq object
its.pd <- merge_phyloseq(asv_its, tax_its, refseq_its, map)
#to_remove_its <- c("DM000013V1PQ83-2", "DM00025V1PQ3141-1", "DM00152V1PQ36-65", "2L75-PE1", "2L66-PE1") #no rpoc data
#its.dat <- prune_samples(!(sample_names(its.pd) %in% to_remove_its), its.pd)
its.dat <- prune_samples(sample_sums(its.pd) > 2000, its.pd) #remove less than 2000 reads
its.dat <- subset_samples(its.dat , hiv_status == "HUU" | hiv_status == "NAN")
its.dat <- subset_taxa(its.dat, V2!="Eukaryote_unknown" & V3!="Fungi_unknown") #remove unwanted taxa

#rpoc
#frequency
seqtab_rpoc <- read.table("~/usa_nigeria/rpoc/sequence_table_rpoc.merged.txt", header=T, row.names=1)
asv_rpoc <-otu_table(seqtab_rpoc, taxa_are_rows=F)
#taxonomy
tax_tab_rpoc <- read.table("~/usa_nigeria/rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
tax_rpoc <- tax_table(as.matrix(tax_tab_rpoc))
#ref sequences
refseq_rpoc <- Biostrings::readDNAStringSet("~/usa_nigeria/rpoc/rep_set_rpoc.fa")
#metadata
metadata <- read.table("~/usa_nigeria/Nigeria_USA_meta.txt", sep="\t", header=T, row.names=1)
map <- sample_data(metadata)
#merge into one phyloseq object
rpoc.pd <- merge_phyloseq(asv_rpoc, tax_rpoc, refseq_rpoc, map)
#to_remove_rpoc <- c("DM00016V1PQ26, DM00078V1PQ32") #not in its data
its_samples <- sample_names(its.dat)
#rpoc.dat <- prune_samples(!(sample_names(rpoc.pd) %in% to_remove_rpoc), rpoc.pd)
rpoc.dat <- prune_samples(sample_sums(rpoc.pd) > 4000, rpoc.pd)
rpoc.dat <- subset_samples(rpoc.dat , hiv_status == "HUU" | hiv_status == "NAN")
rpoc.dat <- prune_samples((sample_names(rpoc.dat) %in% its_samples), rpoc.dat) # get only samples that are in its data as well

rpoc_samples <- sample_names(rpoc.dat)
its.dat <- prune_samples((sample_names(its.dat) %in% rpoc_samples), its.dat) # get only samples that are in its data as well
#merge both phyloseq object 
ps.data <- merge_phyloseq(rpoc.dat, its.dat)

save.image("~/usa_nigeria/phyloseq_obj/ps.RData")