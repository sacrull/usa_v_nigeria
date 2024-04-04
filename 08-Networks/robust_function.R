robustness <- function(g, type=c('vertex', 'edge'),
                       measure=c('btwn.cent', 'degree', 'random'), N=1e3) {
  i <- NULL
  stopifnot(is_igraph(g))
  type <- match.arg(type)
  measure <- match.arg(measure)
  orig_max <-max(components(g)$csize) 
  n <- switch(type, vertex=vcount(g), edge=ecount(g))
  removed.pct <- seq.int(0, 1, length.out=n+1L)
  if (measure == 'random') {
    otype <- paste('Random', type, 'removal')
    rand <- matrix(rep.int(seq_len(n), N), nrow=n, ncol=N)
    index <- apply(rand, 2L, sample)
  } else {
    otype <- paste('Targeted', type, 'attack')
    max.comp.removed <- rep.int(orig_max, n)
  }
  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  if (type == 'vertex') {
    if (measure == 'random') {
      max.comp <- foreach(i=seq_len(N), .combine='cbind') %dopar% {
        ord <- igraph::V(g)$name[index[, i]]
        tmp <- rep.int(orig_max, n)
        g.new <- g
        for (j in seq_len(n - 1L)) {
          g.new <- igraph::delete_vertices(g.new, ord[j])
          tmp[j + 1L] <- max(igraph::components(g.new)$csize)
        }
        tmp
      }
      max.comp.removed <- rowMeans(max.comp)

    } else {
      val <- if (measure == 'btwn.cent') centr_betw(g)$res else check_degree(g)
      ord <- V(g)$name[order(val, decreasing=TRUE)]
      for (j in seq_len(n - 1L)) {
        g <- delete_vertices(g, ord[j])
        max.comp.removed[j + 1L] <- max(components(g)$csize)
      }
    }

  } else {
    if (measure == 'degree') {
      stop('For edge attacks, must choose "btwn.cent" or "random"!')
    } else if (measure == 'random') {
      max.comp <- foreach(i=seq_len(N), .combine='cbind') %dopar% {
        el <- as_edgelist(g, names=FALSE)[index[, i], ]
        tmp <- rep.int(orig_max, n)
        for (j in seq_len(n - 1L)) {
          g.rand <- graph_from_edgelist(el[-seq_len(j), , drop=FALSE], directed=FALSE)
          tmp[j + 1L] <- max(components(g.rand)$csize)
        }
        tmp
      }
      max.comp.removed <- rowMeans(max.comp)

    } else {
      ord <- order(E(g)$btwn, decreasing=TRUE)
      el <- as_edgelist(g, names=FALSE)[ord, ]
      for (j in seq_len(n - 1L)) {
        g <- graph_from_edgelist(el[-seq_len(j), , drop=FALSE], directed=FALSE)
        max.comp.removed[j + 1L] <- max(components(g)$csize)
      }
    }

  }
  max.comp.removed <- c(max.comp.removed, 0)
  comp.pct <- max.comp.removed / orig_max
  out <- data.table(type=otype, measure=measure, comp.size=max.comp.removed,
                    comp.pct=comp.pct, removed.pct=removed.pct)
  if ('name' %in% graph_attr_names(g)) out[, eval(getOption('bg.group')) := g$name]
  return(out)
}