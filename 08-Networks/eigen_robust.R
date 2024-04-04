robustness_eigen <- function(g, type=c('vertex'),
                       measure=c('eigen'), N=1e3) {
i <- NULL
stopifnot(is_igraph(g))
type <- match.arg(type)
measure <- match.arg(measure)
orig_max <-max(components(g)$csize) 
n <- switch(type, vertex=vcount(g), edge=ecount(g))
removed.pct <- seq.int(0, 1, length.out=n+1L)
otype <- paste('Targeted', type, 'attack')
max.comp.removed <- rep.int(orig_max, n)
val <- if (measure == 'eigen') centr_eigen(g)$vector else check_degree(g)
removed.pct <- seq.int(0, 1, length.out=n+1L)

ord <- V(g)$name[order(val, decreasing=TRUE)]
  for (j in seq_len(n - 1L)) {
    g <- delete_vertices(g, ord[j])
    max.comp.removed[j + 1L] <- max(components(g)$csize)
  }
max.comp.removed <- c(max.comp.removed, 0)
comp.pct <- max.comp.removed / orig_max
out <- data.table(type=otype, measure=measure, comp.size=max.comp.removed,
              comp.pct=comp.pct, removed.pct=removed.pct)
if ('name' %in% graph_attr_names(g)) out[, eval(getOption('bg.group')) := g$name]
return(out)
}