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
install.packages("tidyquant")
install.packages("ggdist")
install.packages("ggthemes")
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
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr")+
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
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr")+
  geom_jitter(aes(color=gender), shape=16, position=position_jitter(0.2), size=2.5)
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
