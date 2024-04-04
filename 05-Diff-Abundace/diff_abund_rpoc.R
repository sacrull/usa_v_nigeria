library(phyloseq)
library(corncob)
library(magrittr)

set.seed(12349)
setwd("~/usa_nigeria/diff_abund")
load("~/usa_nigeria/phyloseq_obj/ps.RData")
#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8")# %>% subset_samples(Tooth_Classification == "H-CF") 
rpoc.dat.sp 
rpoc.dat.sub <- rpoc.dat %>% subset_samples(Tooth_Classification == "H-CF") 

#choose a model
corncob <- bbdml(formula = ASV1 ~ 1,
phi.formula = ~ 1,
data = rpoc.dat.sp)

corncob_da <- bbdml(formula = ASV1 ~ Geog_loc,
phi.formula = ~ Geog_loc,
data = rpoc.dat.sp)

lrtest(mod_null = corncob, mod = corncob_da) #got a p-value of less than 0.05 -> want to use covariate model

#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc

#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.geo_loc.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
# differential variance test
dv_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ 1,
                               test = "LRT",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
dv_analysis_rpoc$significant_taxa
pdf("./diffvar.geo_loc.rpoc.pdf", width = 20)
plot(dv_analysis_rpoc, level=c("V8"))
dev.off()

#HCF
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Tooth_Classification == "H-CF")
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc

#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.geo_loc.hcf.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#H
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(tooth_health == "H")
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc

#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.geo_loc.h.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#CF
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(mouth_health == "CF")
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc

#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.geo_loc.cf.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#DCD
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Tooth_Classification == "D-CD") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
da_analysis_rpoc$significant_taxa
pdf("./diffab.geo_loc.dcd.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#CD
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(mouth_health == "CD") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
da_analysis_rpoc$significant_taxa
pdf("./diffab.geo_loc.cd.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#D
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(tooth_health == "D") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
da_analysis_rpoc$significant_taxa
pdf("./diffab.geo_loc.d.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()


#health
#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Tooth_Classification == "H-CF" | Tooth_Classification == "D-CD") 
da_analysis_rpoc <- differentialTest(formula = ~ Tooth_Classification,
                               phi.formula = ~ Tooth_Classification +Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Tooth_Classification +Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)        
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.tooth_class.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#mouth health
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(mouth_health == "CF" | mouth_health == "CD") 
da_analysis_rpoc <- differentialTest(formula = ~ mouth_health,
                               phi.formula = ~ mouth_health+Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ mouth_health+Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)        
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.mouth.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#tooth health
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(tooth_health == "H" | tooth_health == "D") 
da_analysis_rpoc <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health+Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health+Geog_loc,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)        
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.tooth_health.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()
#health categories by geo location
#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "USA") %>% subset_samples(Tooth_Classification == "H-CF" | Tooth_Classification == "D-CD") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ Tooth_Classification,
                               phi.formula = ~ Tooth_Classification,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Tooth_Classification,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.tooth.usa.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()

#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "Nigeria") %>% subset_samples(Tooth_Classification == "H-CF" | Tooth_Classification == "D-CD") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ Tooth_Classification,
                               phi.formula = ~ Tooth_Classification,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Tooth_Classification,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.tooth.nigeria.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()

#mouth health
#health categories by geo location
#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "USA") %>% subset_samples(mouth_health == "CF" | mouth_health == "CD") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ mouth_health,
                               phi.formula = ~ mouth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ mouth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.mouth.usa.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()

#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "Nigeria") %>% subset_samples(mouth_health == "CF" | mouth_health == "CD") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ mouth_health,
                               phi.formula = ~ mouth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ mouth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.mouth.nigeria.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()

#tooth health
#health categories by geo location
#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "USA") %>% subset_samples(tooth_health == "H" | tooth_health == "D") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.tooth_health.usa.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()

#glom to species level
rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Geog_loc == "Nigeria") %>% subset_samples(tooth_health == "H" | tooth_health == "D") 
rpoc.dat.sp 
#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ tooth_health,
                               phi.formula = ~ tooth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ tooth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = rpoc.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_rpoc
#look at sign taxa
da_analysis_rpoc$significant_taxa
pdf("./diffab.tooth_health.nigeria.rpoc.pdf", width = 20)
plot(da_analysis_rpoc, level=c("V8"))
dev.off()

