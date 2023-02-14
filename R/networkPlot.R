#'
#' Enrichment Network
#'
#' @description Displays enrichment results as an easy to interpret cluster network.
#'
#' @param enrichment a data frame containing enrichment results.
#' @param simMethod method for calculating similarity matrix between the pathways. Available
#' methods: \code{'jaccard'}, \code{'cosine'} and \code{'cor'}
#' @param clustMethod method for detecting pathway clusters. Available methods: \code{'markov'},
#' \code{'hier'} and \code{'spectral'}. Using \code{'spectral'} method requires that you have
#' \code{Spectrum} package installed.
#' @param clustNameMethod method for selecting cluster names. Available methods: \code{'pagerank'},
#' \code{'hits'} and \code{'none'}
#' @param colorBy which column should be used to color the nodes. If NULL, will be selected by
#' default (either NES or p-value will be used)
#' @param nodeSize which column should be used to set the node size. If NULL, will be selected by
#' default (either setSize or Count will be used)
#' @param innerCutoff similarity cutoff for in-cluster nodes
#' @param outerCutoff similarity cutoff for between-cluster nodes
#' @param colorType how to colour the nodes: \code{'nes'} - will center around 0 with blue min and
#' red max, \code{'pval'} - will use log transform on the colorBy column and adjust color range
#' @param pCutoff adjust p-value colouring cutoff when using \code{colorType = 'pval'}
#' @param drawEllipses enable / disable ellipse drawing
#' @param fontSize adjust cluster label font size
#' @param repelLabels should cluster label positions be corrected
#' @param minClusterSize min number of nodes in a single cluster
#' @param plotOnly default: TRUE. If set to false, returns clustering results as well as the plot
#' (a list object is returned).
#' @param verbose enable / disable log messages
#'
#' @return \code{ggplot} object.
#'
#' @seealso \code{enrichmentData}, \code{validateEnrichment}
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
  innerCutoff = 0.1,
  outerCutoff = 0.5,
  colorType = c('nes', 'pval'),
  pCutoff = -10,
  drawEllipses = FALSE,
  fontSize = 5,
  repelLabels = FALSE,
  minClusterSize = 2,
  plotOnly = TRUE,
  verbose = FALSE
) {
  if (class(enrichment) != 'data.frame') {
    stop('An object of class data.frame is expected.')
  }

  simMethod <- match.arg(simMethod)
  clustMethod <- match.arg(clustMethod)
  clustNameMethod <- match.arg(clustNameMethod)
  colorType <- match.arg(colorType)

  # Validate the parameters
  if (!between(innerCutoff, 0, 1)) {
    stop('innerCutoff must be between 0 and 1.')
  }

  if (!between(outerCutoff, 0, 1)) {
    stop('outerCutoff must be between 0 and 1.')
  }

  if (minClusterSize < 2) {
    stop('Currently supported minClusterSize >= 2')
  }

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
                           minClusterSize = minClusterSize,
                           verbose = verbose)

  enrichClust <- enrichmentNetwork.prepareEnrichmentClusters(enrichment, clusters, params)

  plot <- enrichmentNetwork.plot(enrichClust,
                                 sim,
                                 clusters,
                                 innerCutoff = innerCutoff,
                                 outerCutoff = outerCutoff,
                                 colorType = colorType,
                                 pCutoff = pCutoff,
                                 drawEllipses = drawEllipses,
                                 fontSize = fontSize,
                                 repelLabels = repelLabels)

  if (plotOnly) {
    return(plot)
  } else {
    return(list(plot = plot,
                clusters = enrichClust[ , list(ID, Cluster) ]))
  }
}

#'
#' Enrichment Clusters
#'
#' @description Prepares enrichment and clustering results for \code{enrichmentNetwork.plot} method.
#'
#' @import data.table
#' @importFrom dplyr %>%
#'
enrichmentNetwork.prepareEnrichmentClusters <- function(enrichment, clusters, params) {
  clusterSizes <- table(clusters)

  data.table(
    ID = enrichment[ , 'Description' ],
    size = enrichment[ , params$nodeSize ],
    color = enrichment[ , params$colorBy ]
  ) %>%
    merge(data.table(Cluster = clusters, ID = names(clusters)), by.x = 'ID', by.y = 'ID') %>%
    .[ , `Cluster size` := as.integer(clusterSizes[ Cluster ]) ]
}

#'
#' Enrichment Network
#'
#' Creates a \code{ggplot} object from pathway similarity matrix and clustering results.
#'
#' @param dt output from \code{enrichmentNetwork.prepareEnrichmentClusters}
#' @param sim similarity matrix
#' @param clust assigned clusters
#' @param innerCutoff similarity cutoff for in-cluster nodes
#' @param outerCutoff similarity cutoff for between-cluster nodes
#' @param colorType how to colour the nodes: \code{'nes'} - will center around 0 with blue min and
#' red max, \code{'pval'} - will use log transform on the colorBy column and colour range [-10, 0]
#' @param pCutoff adjust p-value colouring cutoff when using \code{colorType = 'pval'}
#' @param drawEllipses enable / disable ellipse drawing
#' @param fontSize adjust cluster label font size
#' @param repelLabels should cluster label positions be corrected
#'
#' @return \code{ggplot} object.
#'
#' @import igraph
#' @import data.table
#' @import ggplot2
#' @import ggforce
#'
enrichmentNetwork.plot <- function(dt,
                                   sim,
                                   clust,
                                   innerCutoff = 0.1,
                                   outerCutoff = 0.5,
                                   colorType = c('nes', 'pval'),
                                   pCutoff = -10,
                                   drawEllipses = FALSE,
                                   fontSize = 5,
                                   repelLabels = FALSE
) {
  colorType <- match.arg(colorType)

  graph <- enrichmentNetwork.connect(sim, clust, innerCutoff = innerCutoff, outerCutoff = outerCutoff)
  coordinates <- merge(graph$coordinates, dt, by.x = 'ID', by.y = 'ID')

  lines <- graph$edges
  lines <- lines[ from %in% coordinates[ , ID ] & to %in% coordinates[ , ID ] ]

  if (colorType == 'nes') {
    range <- max(abs(coordinates[ , color ]))
    colors <- scale_color_distiller(limits = c(-range, range), palette = 'Spectral')
    colorTitle <- 'NES'
  }

  if (colorType == 'pval') {
    coordinates[ , color := log(color) ] %>%
      .[ color < pCutoff, color := pCutoff ]
    colors <- scale_color_distiller(limits = c(pCutoff, 0), direction = -1, palette = 'OrRd')
    colorTitle <- 'p-value'
  }

  plot <- ggplot()

  if (drawEllipses) {
    plot <- enrichmentNetwork.addEllipses(plot, coordinates)
  }

  plot <- plot +
    geom_link0(data = lines, aes(x = xStart, y = yStart, xend = xEnd, yend = yEnd), size = 0.25, alpha = 0.2) +
    geom_point(data = coordinates, aes(x = x, y = y, ID = ID, color = color, size = size, Cluster = Cluster, `Cluster size` = `Cluster size`)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = 'right',
          legend.key = element_rect(fill = 'white')) +
    labs(color = colorTitle, size = 'Pathway size') +
    coord_fixed() +
    enrichmentNetwork.clusterLabels(coordinates, fontSize, repelLabels) +
    colors

  plot
}

#'
#' Enrichment Network
#'
#' @description Connects related pathways.
#'
#' @param sim similarity matrix
#' @param clusters a list of clusters
#' @param innerCutoff similarity cutoff for in-cluster nodes
#' @param outerCutoff similarity cutoff for between-cluster nodes
#'
#' @return A list with two objects: \code{coordinates} - each pathway coordinates in the graph, and
#' \code{edges} - list of edges in the graph, with their coordinates.
#'
#' @import igraph
#'
enrichmentNetwork.connect <- function(sim, clusters, innerCutoff = 0.1, outerCutoff = 0.5) {
  stopifnot(rownames(sim) == colnames(sim))
  paths <- rownames(sim)

  counter <- 1
  edges <- list()

  for (i in 1:(nrow(sim) - 1)) {
    for (j in (i + 1):ncol(sim)) {
      value <- sim[ i, j ]

      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        if (value >= innerCutoff) {
          edges[[counter]] <- data.table(from = paths[ i ], to = paths[ j ], weight = value)
          counter <- counter + 1
        }
      }

      if (value >= outerCutoff) {
        edges[[counter]] <- data.table(from = paths[ i ], to = paths[ j ], weight = value)
        counter <- counter + 1
      }
    }
  }

  edges <- rbindlist(edges)

  g <- graph_from_data_frame(edges, directed = FALSE)

  coordinates <- layout_nicely(g) %>% as.data.table
  colnames(coordinates) <- c('x', 'y')
  coordinates[ , ID := V(g) %>% as.list %>% names ]

  edges %>%
    .[ , xStart := lapply(from, \(path) coordinates[ ID %in% path, x ]) %>% unlist ] %>%
    .[ , yStart := lapply(from, \(path) coordinates[ ID %in% path, y ]) %>% unlist ] %>%
    .[ , xEnd := lapply(to, \(path) coordinates[ ID %in% path, x ]) %>% unlist ] %>%
    .[ , yEnd := lapply(to, \(path) coordinates[ ID %in% path, y ]) %>% unlist ]

  list(coordinates = coordinates,
       edges = edges)
}

#'
#' Enrichment Network
#'
#' @description Adds ellipses around pathway clusters.
#'
#' @param plot \code{ggplot} object
#' @param pathways a data table containing pathway coordinates in the graph and their clustering
#' results
#'
#' @import data.table
#'
enrichmentNetwork.addEllipses <- function(plot, pathways) {
  for (cluster in pathways[ , unique(Cluster) ]) {
    points <- pathways[ Cluster == cluster, list(x, y) ]

    if (nrow(points) == 1) { next }
    if (nrow(points) == 2) {
      points <- rbind(points, data.table(x = mean(points[ , x ]) + 0.1, y = mean(points[ , y ]) - 0.1))
    }

    params <- mvee(points, plotme = FALSE)

    ellipsePoints <- enrichmentNetwork.ellipse(params$c, params$a, params$b, params$alpha, bigger = 0.35)

    plot <- plot + geom_path(data = ellipsePoints, aes(x, y), alpha = 0.4, size = 0.2)
  }

  plot
}

#'
#' Enrichment Network
#'
#' @description Calculates ellipse coordinates from its equation.
#'
#' @param c ellipse center
#' @param a a axis
#' @param b b axis
#' @param alpha ellipse angle
#' @param bigger how much to increase ellipse size
#'
#' @import data.table
#'
enrichmentNetwork.ellipse <- function(c, a, b, alpha, bigger = 0) {
  theta <- seq(0, 2 * pi, length = 101)
  a <- a + bigger
  b <- b + bigger

  x <- c[ 1 ] + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
  y <- c[ 2 ] +
    a * cos(theta) * sin(alpha) +
    b * sin(theta) * cos(alpha)

  data.table(x = x, y = y)
}

#'
#' Cluster Labels
#'
#' @description Displays a label for each pathway cluster in the graph.
#'
#' @param pathways a data table containing pathway coordinates in the graph and their clustering
#' results
#' @param fontSize adjust cluster label font size
#' @param repelLabels should cluster label positions be corrected
#'
#' @return A \code{geom_text} object with labels and their coordinates, ready to add to a
#' \code{ggplot} object.
#'
#' @import data.table
#' @import foreach
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr %>%
#'
enrichmentNetwork.clusterLabels <- function(pathways, fontSize = 5, repelLabels = FALSE) {
  if (repelLabels) {
    labels <- foreach(cluster = pathways[ , unique(Cluster) ], .combine = rbind) %do% {
      points <- pathways[ Cluster == cluster, list(x, y) ]
      midPoint <- list(
        x = points[ , x ] %>% mean,
        y = points[ , y ] %>% max
      )

      data.table(x = midPoint$x, y = midPoint$y, label = splitWords(cluster))
    }

    geom_text_repel(data = labels, aes(x = x, y = y, label = label), size = fontSize, segment.size = 0.25)
  } else {
    labels <- foreach(cluster = pathways[ , unique(Cluster) ], .combine = rbind) %do% {
      points <- pathways[ Cluster == cluster, list(x, y) ]
      midPoint <- list(
        x = points[ , x ] %>% mean,
        y = points[ , y ] %>% max + 0.25
      )

      data.table(x = midPoint$x, y = midPoint$y, label = splitWords(cluster))
    }

    geom_text(data = labels, aes(x = x, y = y, label = label), size = fontSize, segment.size = 0.25)
  }
}
