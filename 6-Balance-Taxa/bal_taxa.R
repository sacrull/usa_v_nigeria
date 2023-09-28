#devtools::install_github(repo = "malucalle/selbal")

library(selbal)
library(grid)

setwd("~/usa_nigeria/bal_tax")
load("~/usa_nigeria/phyloseq_obj/ps.RData")

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(its.dat, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.10*length(x)), TRUE)
glom <- microbiome::transform(glom, 'compositional')
# pull data
dat <- as.data.frame(otu_table(glom))
map <- sample_data(glom)

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,1:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# response variable
dif2 <- dim(dat)[2] - 1
y <- factor(dat[,dif2]) #geog lcaotion
z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound


plz <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 10,
                           covar = NULL, logit.acc = "AUC")
# using tooth type as the covariate
bal <- selbal(x = x, y = y, n.fold = 5, n.iter = 10,
                          covar = NULL)
selbal.cv(x = x, y = y, n.fold = 5, n.iter = 10,
                           covar = NULL, logit.acc = "AUC")
# run selbal.cv as well
cv.bal <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 100, covar = z, zero.rep = "bayes", seed = 4597)
cv.bal$global.balance
# plot correlation
pdf("global_plot.pdf")
grid.draw(bal$global.plot)
dev.off()
pdf("global_plot.cv.pdf")
grid.draw(cv.bal$global.plot)
dev.off()
pdf("var_barplot.cv.pdf")
cv.bal$var.barplot
dev.off()
pdf("plot_tab.cv.pdf")
plot.tab(cv.bal$cv.tab)
dev.off()
```

What if we don't control for tooth type?

```R
# using tooth type as the covariate
bal <- selbal(x = x, y = y, zero.rep = "bayes")
# run selbal.cv as well
cv.bal <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 100, zero.rep = "bayes", seed = 4597)
cv.bal$global.balance
# plot correlation
pdf("global_plot.nocovar.pdf")
grid.draw(bal$global.plot)
dev.off()
pdf("global_plot.cv.nocovar.pdf")
grid.draw(cv.bal$global.plot)
dev.off()
pdf("var_barplot.cv.nocovar.pdf")
cv.bal$var.barplot
dev.off()
pdf("plot_tab.cv.nocovar.pdf")
plot.tab(cv.bal$cv.tab)
dev.off()