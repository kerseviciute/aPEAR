#'
#' Enrichment Network
#'
#' @description Displays enrichment results as an easy to interpret cluster network.
#'
#' @param enrichment a data frame containing enrichment results.
#' @param simMethod method for calculating similarity matrix between the pathways. Available
#' methods: \code{'jaccard'}, \code{'cosine'} and \code{'cor'}
#' @param clustMethod method for detecting pathway clusters. Available methods: \code{'markov'},
#' \code{'hier'} and \code{'spectral'}
#' @param clustNameMethod method for selecting cluster names. Available methods: \code{'pagerank'},
#' \code{'hits'} and \code{'none'}
#' @param colorBy which column should be used to color the nodes. If NULL, will be selected by
#' default (either NES or p-value will be used)
#' @param nodeSize which column should be used to set the node size. If NULL, will be selected by
#' default (either setSize or Count will be used)
#' @param verbose enable / disable log messages
#'
#' @seealso \code{enrichmentData}
#'
#' @export
#'
enrichmentNetwork <- function(
  enrichment,
  simMethod = c('jaccard', 'cosine', 'cor'),
  clustMethod = c('markov', 'hier', 'spectral'),
  clustNameMethod = c('pagerank', 'hits', 'none'),
  colorBy = NULL,
  nodeSize = NULL,
  verbose = FALSE
) {
  if (class(enrichment) != 'data.frame') {
    stop('An object of class data.frame is expected.')
  }

  simMethod <- match.arg(simMethod)
  clustMethod <- match.arg(clustMethod)
  clustNameMethod <- match.arg(clustNameMethod)

  params <- validateEnrichment(enrichment,
                               colorBy = colorBy,
                               nodeSize = nodeSize,
                               verbose = verbose)

  sim <- pathwaySimilarity(enrichment,
                           geneCol = params$genesCol,
                           method = simMethod)

  clusters <- findClusters(sim,
                           method = clustMethod,
                           nameMethod = clustNameMethod,
                           verbose = verbose)
  #
  # p <- drawPlot(enrichment, sim, clusters, colorBy, nodeSize)
  #
  # p
}

#'
#' @import data.table
#' @importFrom dplyr %>%
#' @import ggplot2
#' @import ggforce
#' @import reticulate
#' @import foreach
drawPlot <- function(enrichment, sim, clusters, colorBy, innerCutoff = 0.1, outerCutoff = 0.5) {
  # TODO: THIS IS TERRIBLE! FIND A WAY TO MAKE IT WORK WITHOUT USING PYTHON!
  use_condaenv('snakemake')
  nx <- import("networkx")

  graph <- connect(sim, clusters, innerCutoff = innerCutoff, outerCutoff = outerCutoff)
  pos <- nx$spring_layout(graph, seed = as.integer(1234))

  clusterSizes <- table(clusters)

  res <- foreach(id = names(pos), .combine = rbind) %do% {
    data.table(ID = id, x = pos[[id]][ 1 ], y = pos[[id]][ 2 ])
  } %>%
    merge(enrichment[ , c('Description', colorBy, 'setSize') ], by.x = 'ID', by.y = 'Description') %>%
    merge(data.table(Cluster = clusters, ID = names(clusters)), by.x = 'ID', by.y = 'ID') %>%
    .[ , `Cluster size` := as.integer(clusterSizes[ Cluster ]) ]

  range <- max(abs(res[ , NES ]))
  lines <- listEdges(graph, res)

  plot <- ggplot()

  plot <- addEllipses(plot, res, clusters)

  plot <- plot +
    geom_link0(data = lines, aes(x = x, y = y, xend = xend, yend = yend), size = 0.1, alpha = 0.3) +
    geom_point(data = res, aes(x = x, y = y, ID = ID, color = NES, size = setSize, Cluster = Cluster, `Cluster size` = `Cluster size`)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    scale_color_distiller(limits = c(-range, range), palette = 'Spectral') +
    coord_fixed() +
    clusterLabels(res) +
    ylim(-1.1, 1.1) +
    xlim(-1.1, 1.1)
}

#' @import reticulate
connect <- function(sim, clusters, innerCutoff = 0.1, outerCutoff = 0.5) {
  use_condaenv('snakemake')
  nx <- import("networkx")

  stopifnot(rownames(sim) == colnames(sim))
  paths <- rownames(sim)
  g <- nx$Graph()

  for (i in 1:(nrow(sim) - 1)) {
    for (j in (i + 1):ncol(sim)) {
      value <- sim[ i, j ]

      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        if (value >= innerCutoff) g$add_edge(paths[ i ], paths[ j ], weight = value)
        next
      }

      if (value >= outerCutoff) {
        g$add_edge(paths[ i ], paths[ j ], weight = value)
      }
    }
  }

  return(g)
}

#' @import reticulate
#' @import data.table
listEdges <- function(graph, positions) {
  iter <- as_iterator(graph$edges)
  current <- iter_next(iter)

  edges <- list()
  counter <- 0
  while (!is.null(current)) {
    x <- current

    path1 <- x[[1]]
    if (!(path1 %in% positions[ , ID ])) {
      current <- iter_next(iter)
      next
    }
    path1 <- positions[ ID == path1 ]

    path2 <- x[[2]]
    if (!(path2 %in% positions[ , ID ])) {
      current <- iter_next(iter)
      next
    }
    path2 <- positions[ ID == path2 ]

    counter <- counter + 1
    edges[[counter]] <- data.table(x = path1[ , x ], y = path1[ , y ], xend = path2[ , x ], yend = path2[ , y ])

    current <- iter_next(iter)
  }

  return(rbindlist(edges))
}

#' @import data.table
addEllipses <- function(plot, pathways, clusters) {
  for (cluster in unique(clusters)) {
    points <- pathways[ Cluster == cluster, list(x, y) ]

    if (nrow(points) == 1) { next }
    if (nrow(points) == 2) {
      points <- rbind(points, data.table(x = mean(points[ , x ]) + 0.01, y = mean(points[ , y ]) - 0.01))
    }

    params <- mvee(points, plotme = FALSE)
    if (params$ab[ 1 ] > 0.5 | params$ab[ 2 ] > 0.5) next

    ellipsePoints <- ellipse(params$c, params$ab, params$alpha, bigger = 0.05)

    plot <- plot + geom_path(data = ellipsePoints, aes(x, y), alpha = 0.4, size = 0.2)
  }

  plot
}

ellipse <- function(c, ab, alpha, bigger = 0) {
  pacman::p_load(data.table)
  theta <- seq(0, 2 * pi, length = 101)
  a <- ab[ 1 ] + bigger
  b <- ab[ 2 ] + bigger

  x <- c[ 1 ] + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
  y <- c[ 2 ] +
    a * cos(theta) * sin(alpha) +
    b * sin(theta) * cos(alpha)

  data.table(x = x, y = y)
}

#' @import data.table
clusterLabels <- function(pathways) {
  library(foreach)
  library(data.table)
  library(dplyr)
  library(ggplot2)

  labels <- foreach(cluster = pathways[ , unique(Cluster) ], .combine = rbind) %do% {
    points <- pathways[ Cluster == cluster, list(x, y) ]
    midPoint <- list(
      x = points[ , x ] %>% mean,
      y = points[ , y ] %>% mean + 0.1
    )

    data.table(x = midPoint$x, y = midPoint$y, label = splitWords(cluster))
  }

  geom_text(data = labels, aes(x = x, y = y, label = label), size = 2)
}

splitWords <- function(x) {
  strwrap(x, width = 30) %>% paste(collapse = '\n')
}
