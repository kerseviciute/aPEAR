#' 
#' Pathway Clustering
#' 
#' @description Finds clusters in a similarity matrix.
#' 
#' @param sim similarity matrix
#' @param minClusterSize minimum cluster size (integer, >= 1). This parameter is ignored when using
#' \code{method = 'spectral'}
#' @param method a method to be used for detecting clusters. Available methods are: \code{'markov'},
#' \code{'hier'} and \code{'spectral'}. Using \code{'spectral'} method requires that you have
#' \code{Spectrum} package installed.
#' @param nameMethod a method for setting cluster names. Available methods are: \code{'pagerank'},
#' \code{'hits'}, and \code{'none'}
#' @param verbose enable / disable log messages
#' 
#' @importFrom dplyr %>%
#' 
#' @export
#' 
findClusters <- function(sim,
                         minClusterSize = 2,
                         method = c('markov', 'hier', 'spectral'),
                         nameMethod = c('pagerank', 'hits', 'none'),
                         verbose = FALSE
) {
  method <- match.arg(method)
  nameMethod <- match.arg(nameMethod)

  if (minClusterSize < 1) {
    warning(paste0('Invalid minClusterSize = ', minClusterSize, '. Will set minClusterSize = 2.'))
    minClusterSize <- 2
  }

  clusters <- switch(method,
                     'markov' = findClustersMarkov(sim, minClusterSize, verbose),
                     'hier' = findClustersHier(sim, minClusterSize, verbose),
                     'spectral' = findClustersSpectral(sim, verbose))

  if (length(clusters) == 0) {
    stop('No clusters found.')
  }

  findClusterNames(sim, clusters, nameMethod)
}

#' 
#' Markov Cluster Algorithm
#' 
#' @description Finds clusters in a similarity matrix using Markov Cluster Algorithm.
#' 
#' @param sim similarity matrix
#' @param minClusterSize minimum cluster size
#' @param verbose enable / disable log messages
#' 
#' @importFrom MCL mcl
#' @importFrom dplyr %>%
#' 
findClustersMarkov <- function(sim,
                               minClusterSize = 2,
                               verbose = FALSE
) {
  if (verbose) message('Using Markov Cluster Algorithm to detect pathway clusters...')

  allow1 <- ifelse(minClusterSize == 1, TRUE, FALSE)

  res <- mcl(sim,
             addLoops = FALSE,
             max.iter = 500,
             expansion = 2,
             inflation = 2.5,
             allow1 = allow1)

  if (!('Cluster' %in% names(res))) {
    stop('Unable to cluster data using Markov clustering algorithm.')
  }

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

  if (verbose) message('Clustering done')

  clusters
}

#' 
#' Hierarchical Clustering
#' 
#' @description Finds clusters using hierarchical algorithm.
#' 
#' @param sim similarity matrix
#' @param minClusterSize minimum cluster size
#' @param verbose enable / disable log messages
#' 
findClustersHier <- function(sim,
                             minClusterSize = 2,
                             verbose = FALSE
) {
  if (verbose) message('Using Hierarchical Clustering to detect pathway clusters...')

  hCluster <- as.dist(1 - sim) %>% hclust
  clusters <- cutree(hCluster, h = 0.9)

  if (minClusterSize != 1) {
    sizes <- clusters %>% table
    ids <- names(sizes)[ sizes >= minClusterSize ]
    clusters <- clusters[ clusters %in% ids ]
  }

  if (verbose) message('Clustering done')

  clusters
}

#' 
#' Spectral Clustering
#' 
#' @description Finds clusters using adaptive spectral clustering. Using this method requires that
#' you have \code{Spectrum} package installed.
#' 
#' @param sim similarity matrix
#' @param verbose enable / disable log messages
#' 
findClustersSpectral <- function(sim, verbose = FALSE) {
  if (verbose) message('Using Spectral Clustering to detect pathway clusters...')
  library(Spectrum)

  clusters <- Spectrum(sim, maxk = 100, showres = FALSE, silent = !verbose, clusteralg = 'km')
  clusters <- clusters$assignments
  names(clusters) <- colnames(sim)

  if (verbose) message('Clustering done')

  clusters
}
