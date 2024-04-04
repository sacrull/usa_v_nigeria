library(phyloseq)
library(corncob)
library(magrittr)

set.seed(12349)
setwd("~/usa_nigeria/diff_abund")
load("~/usa_nigeria/phyloseq_obj/ps.RData")
#glom to species level
its.dat <- prune_samples(sample_sums(its.pd) > 2000, its.pd) %>% subset_samples(hiv_status == "HUU" | hiv_status == "NAN") %>%
  subset_taxa(V2!="Eukaryote_unknown" & V3!="Fungi_unknown") #remove unwanted taxa
its.dat.sp <- its.dat %>%
                tax_glom("V8")# %>% subset_samples(Tooth_Classification == "H-CF") 
its.dat.sp 
its.dat.sub <- its.dat %>% subset_samples(Tooth_Classification == "H-CF") 

#choose a model
corncob <- bbdml(formula = ASV1 ~ 1,
phi.formula = ~ 1,
data = its.dat.sp)

corncob_da <- bbdml(formula = ASV1 ~ Geog_loc,
phi.formula = ~ Geog_loc,
data = its.dat.sp)

lrtest(mod_null = corncob, mod = corncob_da) #got a p-value of less than 0.05 -> want to use covariate model

#diff abundance
da_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its

#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.geo_loc.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
# differential variance test
dv_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ 1,
                               test = "LRT",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
dv_analysis_its$significant_taxa
pdf("./diffvar.geo_loc.its.pdf", width = 20)
plot(dv_analysis_its, level=c("V8"))
dev.off()

#HCF
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Tooth_Classification == "H-CF")
da_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its

#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.geo_loc.hcf.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#H
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(tooth_health == "H")
da_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its

#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.geo_loc.h.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#CF
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(mouth_health == "CF")
da_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its

#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.geo_loc.cf.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#DCD
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Tooth_Classification == "D-CD") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
da_analysis_its$significant_taxa
pdf("./diffab.geo_loc.dcd.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#CD
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(mouth_health == "CD") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
da_analysis_its$significant_taxa
pdf("./diffab.geo_loc.cd.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#D
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(tooth_health == "D") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
da_analysis_its$significant_taxa
pdf("./diffab.geo_loc.d.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()


#health
#glom to species level
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Tooth_Classification == "H-CF" | Tooth_Classification == "D-CD") 
da_analysis_its <- differentialTest(formula = ~ Tooth_Classification,
                               phi.formula = ~ Tooth_Classification,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Tooth_Classification,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)        
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.tooth_class.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#mouth health
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(mouth_health == "CF" | mouth_health == "CD") 
da_analysis_its <- differentialTest(formula = ~ mouth_health,
                               phi.formula = ~ mouth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ mouth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)        
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.mouth.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#tooth health
its.dat.sub <- subset_samples(its.dat.sp, tooth_health == "H" | tooth_health == "D") 
da_analysis_its <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sub,
                               fdr_cutoff = 0.05)        
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.tooth_health.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
#health categories by geo location
#glom to species level
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "USA") %>% subset_samples(Tooth_Classification == "H-CF" | Tooth_Classification == "D-CD") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ Tooth_Classification,
                               phi.formula = ~ Tooth_Classification,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Tooth_Classification,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.tooth.usa.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()

#glom to species level
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "Nigeria") %>% subset_samples(Tooth_Classification == "H-CF" | Tooth_Classification == "D-CD") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ Tooth_Classification,
                               phi.formula = ~ Tooth_Classification,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Tooth_Classification,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.tooth.nigeria.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()

#mouth health
#health categories by geo location
#glom to species level
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "USA") %>% subset_samples(mouth_health == "CF" | mouth_health == "CD") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ mouth_health,
                               phi.formula = ~ mouth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ mouth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.mouth.usa.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()

#glom to species level
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "Nigeria") %>% subset_samples(mouth_health == "CF" | mouth_health == "CD") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ mouth_health,
                               phi.formula = ~ mouth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ mouth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.mouth.nigeria.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()

#tooth health
#health categories by geo location
#glom to species level
its.dat.sp <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "USA") %>% subset_samples(tooth_health == "H" | tooth_health == "D") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.tooth_health.usa.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()

#glom to species level
its.dat.sub <- its.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "Nigeria") %>% subset_samples(tooth_health == "H" | tooth_health == "D") 
its.dat.sp 
#diff abundance
da_analysis_its <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sub,
                               fdr_cutoff = 0.05)
da_analysis_its
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.tooth_health.nigeria.its.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()

