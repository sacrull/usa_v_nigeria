library(phyloseq)
library(corncob)
library(magrittr)

set.seed(12349)
setwd("~/usa_nigeria/diff_abund")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

rpoc.dat.sp <- rpoc.dat %>%
                tax_glom("V8") %>% subset_samples(Tooth_Classification == "H-CF") 
rpoc.dat.sp 

#choose a model
corncob <- bbdml(formula = rASV1 ~ 1,
phi.formula = ~ 1,
data = rpoc.dat.sp)

corncob_da <- bbdml(formula = rASV1 ~ Geog_loc,
phi.formula = ~ Geog_loc,
data = rpoc.dat.sp)

lrtest(mod_null = corncob, mod = corncob_da) #got a p-value of less than 0.05 -> want to use covariate model

#diff abundance
da_analysis_rpoc <- differentialTest(formula = ~ Geog_loc,
                               phi.formula = ~ Geog_loc,
                               formula_null = ~ 1,
                               phi.formula_null = ~ 1,
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