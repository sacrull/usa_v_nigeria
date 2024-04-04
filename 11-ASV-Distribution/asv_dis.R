library(microbiome)
library(EnvStats)
setwd("~/usa_nigeria/asv_dis")
#rarefied samples
rare_its <- rarefy_even_depth(otu_table(its.dat), rngseed = TRUE, replace = FALSE)
data_rare_its = data.frame(otu_table(rare_its)) # create a separated file
ps.rar.its <- phyloseq(rare_its, tax_its, map) # create a phyloseq object
#debra pros
ps.rar.candida.its <- subset_taxa(its.dat, V8 == "Debaryomyces_prosopidis")
ps.rar.candida.its <- tax_glom(ps.rar.candida.its, "V8" ) #use this for wanting at species level
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

data4 <- filter(data, Abundance> 1)

pdf("./asv_dis.tooth_health.d_pros.pdf")
ggplot(data4,aes(x=tooth_health,y=Abundance, fill =Tooth_Classification)) + 
  geom_bar(stat='identity') +
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
  stat_n_text(data =data2$sample ,color = "black",text.box = TRUE)+
  facet_wrap(~Geog_loc)
dev.off()
#look at signficane
data4 %>%
  # group_by(Geog_loc) %>%
  wilcox_test(Abundance ~ tooth_health) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
pdf("./asv_dis.tooth_health.sig.d_pros.pdf")
ggplot(data4, aes(x=tooth_health,y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=Tooth_Classification), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = healthCols)+ #color dots by sample
  theme_classic()+
  facet_wrap(~Geog_loc)
dev.off()
#rarefied samples
#debra pros
ps.rar.candida.its <- subset_taxa(its.dat, V8 == "Sterigmatomyces_halophilus")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.samples.s_halo.pdf")
ggplot(data,aes(x=tooth_health,y=Abundance, color=Sample)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  facet_wrap(~Geog_loc)
  #scale_fill_manual(values=healthCols)
dev.off()

#debra pros
ps.rar.candida.its <- subset_taxa(its.dat, V8 == "Aspergillus_intermedius")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.samples.a_inter.pdf")
ggplot(data,aes(x=Sample,y=Abundance)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()

#debra pros
ps.rar.candida.its <- subset_taxa(its.dat, V8 == "Candida_tropicalis")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.samples.c_trop.pdf")
ggplot(data,aes(x=Sample,y=Abundance)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()

#debra pros
ps.rar.candida.its <- subset_taxa(its.dat, V8 == "Candida_albicans")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.samples.c_albicans.pdf")
ggplot(data,aes(x=Sample,y=Abundance)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)
dev.off()


#rarefied samples
#debra pros
ps.rar.candida.its <- subset_taxa(its.dat, V8 == "Rhodotorula_mucilaginosa")
ps.rar.candida.its <- tax_glom(ps.rar.candida.its, "V8" ) #use this for wanting at species level
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.candida.its, n=50)
ps_50 <- subset_taxa(ps.rar.its, rownames(tax_table(ps.rar.its)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

data4 <- filter(data, Abundance> 1)

pdf("./asv_dis.tooth_health.r_muc.pdf")
ggplot(data4,aes(x=tooth_health,y=Abundance, fill =Tooth_Classification)) + 
  geom_bar(stat='identity') +
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))+
  scale_fill_manual(values=healthCols)+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
  stat_n_text(data =data2$sample ,color = "black",text.box = TRUE)+
  facet_wrap(~Geog_loc)
dev.off()
#look at signficane
data4 %>%
  # group_by(Geog_loc) %>%
  wilcox_test(Abundance ~ tooth_health) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
pdf("./asv_dis.tooth_health.sig.r_muc.pdf")
ggplot(data4, aes(x=tooth_health,y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=Tooth_Classification), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = healthCols)+ #color dots by sample
  theme_classic()+
  facet_wrap(~Geog_loc)
dev.off()