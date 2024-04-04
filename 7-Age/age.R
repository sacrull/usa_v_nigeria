library(ggplot2, verbose=F)
library(phyloseq, verbose=F)
library(ape, verbose=F)
library(metagMisc, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(vegan, verbose=F)
library(ranacapa, verbose=F)
library(microbiome, verbose=F)
library(corncob, verbose=F)
library(magrittr, verbose=F)
library(ggpubr, verbose=F)
library(ecole, verbose=F)
library(introdataviz, verbose=F)
#install.packages("tidyquant")
#install.packages("ggdist")
#install.packages("ggthemes")
#devtools::install_github("psyteachr/introdataviz")
#set seed
set.seed(12349)
setwd("~/usa_nigeria/age")
load("~/usa_nigeria/phyloseq_obj/ps.RData")
#commmunity composition
#age
sampledf <- data.frame(sample_data(its.dat))
sampledf2<-na.omit(sampledf)
one_anova <- aov(age ~ Geog_loc, data = sampledf2)
summary(one_anova)
pdf("age_comparison.pdf")
ggplot(sampledf2, aes(x=Geog_loc, y=age)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
                outlier.size=4)+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr")+
  geom_jitter(aes(color=Geog_loc), shape=16, position=position_jitter(0.2), size=2.5)
dev.off()


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                        position = "dodge", trim = TRUE, scale = "area",
                        show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
)

my_clrs_yct <- c("#404040", "#407a8c")

level_order = c("Nigeria", "USA","None")
pdf("cloud_age.pdf")
ggplot(sampledf2, aes(x = Aliquot, y = age, fill = Geog_loc)) +
  geom_flat_violin(aes(fill = Geog_loc),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = Geog_loc, y = age, colour = Geog_loc),position = position_jitter(width = .01), size = 1, shape = 20)+
  geom_boxplot(aes(x = Geog_loc, y = age, fill = Geog_loc),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_x_discrete(limits = level_order)+
  coord_flip()+
  theme_classic()
dev.off()

pdf("gender.age.boxplot.pdf")
ggplot(sampledf2, aes(x=gender, y=age)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
                outlier.size=4)+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr")+
  geom_jitter(aes(color=gender), shape=16, position=position_jitter(0.2), size=2.5)
dev.off()




age_control <- subset_samples(its.dat, age <= "6" & age >= "4")
sampledf <- data.frame(sample_data(age_control))
sampledf2<-na.omit(sampledf)
one_anova <- aov(age ~ Geog_loc, data = sampledf2)
summary(one_anova)
pdf("age_comparison_adj.pdf")
ggplot(sampledf2, aes(x=Geog_loc, y=age)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
                outlier.size=4)+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr")+
  geom_jitter(aes(color=Geog_loc), shape=16, position=position_jitter(0.2), size=2.5)
dev.off()

pdf("cloud_age_adj.pdf")
ggplot(sampledf2, aes(x = Aliquot, y = age, fill = Geog_loc)) +
  geom_flat_violin(aes(fill = Geog_loc),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = Geog_loc, y = age, colour = Geog_loc),position = position_jitter(width = .01), size = 1, shape = 20)+
  geom_boxplot(aes(x = Geog_loc, y = age, fill = Geog_loc),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_x_discrete(limits = level_order)+
  coord_flip()+
  theme_classic()
dev.off()


#its differences from age
#clr transform
its.age <- subset_samples(its.dat, age <= "6" & age >= "4")
its.dat.clr <- microbiome::transform(its.age, transform="clr", target="OTU")
#betadiversity
#color pallette
geo_color <- c("#8213A0", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

#tooth healt+geog location
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~Tooth_Classification + Geog_loc)
pdf("./bdiv_cap.tooth_health.age.its.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Tooth_Classification", shape ="Geog_loc") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean") # check signficane
#just geo location
ordcap <- ordinate(its.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo.age.its.pdf")
plot_ordination(its.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(its.dat.clr), grp=sample_data(its.dat.clr)$Geog_loc, method="euclidean") # check signficane

#rpoc differences from age
#clr transform
rpoc.age <- subset_samples(rpoc.dat, age <= "6" & age >= "4")
rpoc.dat.clr <- microbiome::transform(rpoc.age, transform="clr", target="OTU")
#betadiversity
#color pallette
geo_color <- c("#8213A0", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Tooth_Classification + Geog_loc)
pdf("./bdiv_cap.tooth_health.age.rpoc.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Tooth_Classification", shape ="Geog_loc") + 
    theme_minimal()+
    scale_color_manual(values=healthCols)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean") # check signficane
#just geo location
ordcap <- ordinate(rpoc.dat.clr, "CAP", "euclidean", ~Geog_loc)
pdf("./bdiv_cap.geo.age.rpoc.pdf")
plot_ordination(rpoc.dat.clr, ordcap, "samples", color="Geog_loc") + 
    theme_minimal() + 
    coord_fixed() +
    stat_ellipse(type="t", linetype=2)
dev.off()
permanova_pairwise(otu_table(rpoc.dat.clr), grp=sample_data(rpoc.dat.clr)$Geog_loc, method="euclidean") # check signficane











rain_height <- .1

pdf("test.pdf")
ggplot(sampledf2, aes(x = "", y = age, fill = Geog_loc)) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,
    position = position_nudge(x = rain_height+.05)) +
  # rain
  geom_point(aes(colour = Geog_loc), size = 2, alpha = .5, show.legend = FALSE, 
              position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, 
               outlier.shape = NA,
               position = position_nudge(x = -rain_height*2)) +
  # mean and SE point in the cloud
  stat_summary(fun.data = mean_cl_normal, mapping = aes(color = Geog_loc), show.legend = FALSE,
               position = position_nudge(x = rain_height * 3)) +
  # adjust layout
  scale_x_discrete(name = "", expand = c(rain_height*3, 0, 0, 0.7)) +
  scale_y_continuous(name = "Age")+
  coord_flip() +
  # custom colours and theme
  scale_fill_brewer(palette = "Dark2", name = "Language group") +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = "white", color = "white"))
dev.off()

pdf("test.pdf")
ggplot(sampledf2,mapping=aes(x=as.factor(gender)))+
geom_bar(aes(color=gender),width=0.5)
dev.off()

pdf("test.pdf")
ggplot(sampledf2,mapping=aes(x=as.factor(gender)))+
geom_bar(aes(fill=gender),width=0.5)+
facet_wrap(~Geog_loc)+
theme_classic()
dev.off()

pdf("test.pdf")
ggplot(sampledf2,mapping=aes(x=as.factor(Geog_loc)))+
geom_bar(aes(fill=gender),width=0.5)+
facet_wrap(~gender)+
theme_classic()
dev.off()
