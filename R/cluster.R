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
  if (!(nameMethod %in% nameMethods)) {
    stop(paste0('Unrecognized nameMethod "', method, '"'))
  }

  if (minClusterSize < 1) {
    stop(paste0('minClusterSize = ', minClusterSize, '. Please choose minClusterSize >= 1.'))
  }

  if (method == 'markov') return(findClustersMarkov(sim, minClusterSize, nameMethod))
  if (method == 'hier') return(findClustersHier(sim, minClusterSize, nameMethod))
  if (method == 'spectral') return(findClustersSpectral(sim, nameMethod))
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

  findClusterNames(sim, clusters, method)
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

  findClusterNames(sim, clusters, method)
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

  findClusterNames(sim, clusters, method)
}
