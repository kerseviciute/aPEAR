#' Finds clusters in a similarity matrix.
#'
#' @param sim similarity matrix
#' @param minClusterSize minimum cluster size (integer, >= 1). This parameter is ignored when using
#' \code{method = 'spectral'}
#' @param method a method to be used for detecting clusters. Available options include 'markov',
#' 'hier' and 'spectral'
#' @param nameMethod a method for setting cluster names. Available method include 'pagerank',
#' 'hits', 'wordcloud' and 'none'
#'
#' @importFrom dplyr %>%
#'
#' @export
findClusters <- function(sim, minClusterSize = 2, method = 'markov', nameMethod = 'pagerank') {
  methods <- list('markov', 'hier', 'spectral')
  if (!(method %in% methods)) {
    stop(paste0('Unrecognized method "', method, '"'))
  }

  nameMethods <- c('none', 'pagerank', 'wordcloud', 'hits')
  if (!(nameMethods %in% nameMethods)) {
    stop(paste0('Unrecognized nameMethod "', method, '"'))
  }

  if (minClusterSize < 1) {
    stop(paste0('minClusterSize = ', minClusterSize, '. Please choose minClusterSize >= 1.'))
  }

  if (method == 'markov') findClustersMarkov(sim, minClusterSize, nameMethod) %>% return
  if (method == 'hier') findClustersHier(sim, minClusterSize, nameMethod) %>% return
  if (method == 'spectral') findClustersSpectral(sim, nameMethod) %>% return
}

#' Finds clusters in a similarity matrix.
#'
#' @param sim similarity matrix
#' @param minClusterSize minimum cluster size
#' @param method a method for setting the cluster names
#'
#' @importFrom MCL mcl
#' @importFrom dplyr %>%
findClustersMarkov <- function(sim, minClusterSize = 2, method = 'pagerank') {
  allow1 <- ifelse(minClusterSize == 1, TRUE, FALSE)

  res <- mcl(sim,
             addLoops = FALSE,
             max.iter = 500,
             expansion = 2,
             inflation = 2.5,
             allow1 = allow1)

  clusters <- res$Cluster
  names(clusters) <- rownames(sim)

  if (!allow1) {
    clusters <- clusters[ clusters != 0 ]
  }

  if (minClusterSize != 1) {
    sizes <- clusters %>% table
    ids <- names(sizes)[ sizes >= minClusterSize ]
    clusters <- clusters[ clusters %in% ids ]
  }

  setClusterNames(sim, clusters, method)
}

#' Finds clusters using hierarchical algorithm.
#'
#' @param sim similarity matrix
#' @param minClusterSize minimum cluster size
#' @param method a method for setting the cluster names
findClustersHier <- function(sim, minClusterSize = 2, method = 'pagerank') {
  hCluster <- as.dist(1 - sim) %>% hclust
  clusters <- cutree(hCluster, h = 0.9)

  if (minClusterSize != 1) {
    sizes <- clusters %>% table
    ids <- names(sizes)[ sizes >= minClusterSize ]
    clusters <- clusters[ clusters %in% ids ]
  }

  setClusterNames(sim, clusters, method)
}

#' Finds clusters using adaptive spectral clustering.
#'
#' @param sim similarity matrix
#' @param method a method for setting the cluster names
#'
#' @importFrom Spectrum Spectrum
findClustersSpectral <- function(sim, method = 'wordcloud') {
  clusters <- Spectrum(sim, maxk = 50, showres = FALSE)
  clusters <- clusters$assignments
  names(clusters) <- colnames(sim)

  setClusterNames(sim, clusters, method)
}

#' Selects best fitting name for a cluster of pathways.
#'
#' @param sim similarity matrix used to detect the clusters
#' @param clusters a vector of clusters, with names indicating the pathway name
#' @param method a method for setting cluster names. Available method include 'pagerank',
#' 'hits', 'wordcloud' and 'none'
setClusterNames <- function(sim, clusters, method = 'pagerank') {
  if (method == 'pagerank') { return(clusterNamesPagerank(sim, clusters)) }
  if (method == 'wordcloud') { return(clusterNamesWordcloud(clusters)) }
  if (method == 'hits') { return(clusterNamesHits(sim, clusters)) }
  if (method == 'none') { return(clusters) }
}

#' Selects a name for each pathway cluster using word cloud method.
#' TODO: how to use python in R packages?
clusterNamesWordcloud <- function(clusters) {
  stop('wordcloud method not implemented yet')
  pacman::p_load(reticulate)
  pacman::p_load(foreach)
  use_condaenv('snakemake')
  te <- import('ThemeExtractor')

  extractor <- te$ThemeExtractor()

  clusterNames <- foreach(cluster = unique(clusters), .combine = rbind) %do% {
    pathways <- names(clusters)[ clusters == cluster ]
    name <- extractor$extract(pathways, maxWords = 4, maxProb = 0.4)

    data.table(ClusterID = cluster, Name = name)
  }

  as.data.table(clusters, keep.rownames = TRUE) %>%
    merge(clusterNames, by.x = 'clusters', by.y = 'ClusterID') %>%
    .[ , list(rn, Name) ] %>%
    deframe
}

#' Selects a name for each pathway cluster using pagerank method.
#'
#' @param sim a similarity matrix used to detect the clusters
#' @param clusters a vector of clusters, with names indicating the pathway name
#'
#' @import igraph
#' @import data.table
#' @import foreach
#' @importFrom dplyr %>%
clusterNamesPagerank <- function(sim, clusters) {
  stopifnot(rownames(sim) == colnames(sim))

  paths <- rownames(sim)
  edges <- list()
  counter <- 1

  for (i in 1:(nrow(sim) - 1)) {
    for (j in (i + 1):ncol(sim)) {
      value <- sim[ i, j ]

      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        edges[[counter]] <- data.table(from = paths[ i ], to = paths[ j ], weight = value)
        counter <- counter + 1
      }
    }
  }

  edges <- rbindlist(edges)
  g <- graph_from_data_frame(edges, directed = FALSE)
  scores <- page.rank(g)$vector

  clusterNames <- foreach(cluster = unique(clusters), .combine = rbind) %do% {
    name <- scores[ names(scores) %in% names(clusters[ clusters == cluster ]) ] %>%
      which.max %>%
      names

    data.table(ClusterID = cluster, Name = name)
  }

  as.data.table(clusters, keep.rownames = TRUE) %>%
    merge(clusterNames, by.x = 'clusters', by.y = 'ClusterID') %>%
    .[ , list(rn, Name) ] %>%
    deframe
}

#' Selects a name for each pathway cluster using hits method.
#'
#' @param sim a similarity matrix used to detect the clusters
#' @param clusters a vector of clusters, with names indicating the pathway name
#'
#' @importFrom arules hits
#' @import foreach
#' @import data.table
#' @importFrom dplyr %>%
clusterNamesHits <- function(sim, clusters) {
  stopifnot(rownames(sim) == colnames(sim))
  adjacency <- matrix(data = 0, nrow = nrow(sim), ncol = ncol(sim))
  rownames(adjacency) <- rownames(sim)
  colnames(adjacency) <- colnames(sim)

  paths <- rownames(sim)

  for (i in 1:(nrow(sim) - 1)) {
    for (j in (i + 1):ncol(sim)) {
      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        adjacency[ i, j ] <- 1
        adjacency[ j, i ] <- 1
      }
    }
  }

  scores <- hits(adjacency, type = 'relative')

  clusterNames <- foreach(cluster = unique(clusters), .combine = rbind) %do% {
    name <- scores[ names(scores) %in% names(clusters[ clusters == cluster ]) ] %>%
      which.max %>%
      names

    data.table(ClusterID = cluster, Name = name)
  }

  as.data.table(clusters, keep.rownames = TRUE) %>%
    merge(clusterNames, by.x = 'clusters', by.y = 'ClusterID') %>%
    .[ , list(rn, Name) ] %>%
    deframe
}
