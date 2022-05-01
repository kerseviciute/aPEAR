#' 
#' PageRank
#' 
#' @description Selects a name for each pathway cluster using PageRank method.
#' 
#' @param sim a similarity matrix used to detect the clusters
#' @param clusters a vector of clusters, with names indicating the pathway name
#' 
#' @import igraph
#' @import data.table
#' @import foreach
#' @importFrom dplyr %>%
#' 
clusterNamesPagerank <- function(sim, clusters) {
  stopifnot(rownames(sim) == colnames(sim))
  stopifnot(nrow(sim) == ncol(sim))

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

  mapClusterNames(scores, clusters)
}

#' 
#' HITS
#' 
#' @description Selects a name for each pathway cluster using HITS method.
#' 
#' @param sim a similarity matrix used to detect the clusters
#' @param clusters a vector of clusters, with names indicating the pathway name
#' 
#' @importFrom arules hits
#' @import foreach
#' @import data.table
#' @importFrom dplyr %>%
#' 
clusterNamesHits <- function(sim, clusters) {
  adjacency <- adjacencyFromSimilarity(sim, clusters)

  scores <- hits(adjacency, type = 'relative')

  mapClusterNames(scores, clusters)
}

#' 
#' Cluster Name
#' 
#' @description Selects the cluster for each pathway from a precalculated cluster title scores.
#' 
#' @param scores pathway evaluation score as cluster center. It is a list with values as scores
#' and names as pathway descriptions
#' @param clusters a list of clusters, where values are cluster ID and names are pathway names
#' 
#' @import foreach
#' @import data.table
#' @importFrom tibble deframe
#' 
mapClusterNames <- function(scores, clusters) {
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

#' 
#' Adjacency Matrix
#' 
#' @description Creates an adjacency matrix from similarity matrix, connecting nodes that belong to
#' the same cluster.
#' 
#' @param sim similarity matrix
#' @param clusters list of clusters
#' 
adjacencyFromSimilarity <- function(sim, clusters) {
  stopifnot(rownames(sim) == colnames(sim))
  stopifnot(nrow(sim) == ncol(sim))

  d <- nrow(sim)
  paths <- rownames(sim)
  names(paths) <- paths

  adjacency <- emptyMatrix(paths, data = 0)

  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        adjacency[ i, j ] <- 1
        adjacency[ j, i ] <- 1
      }
    }
  }

  adjacency
}

#' 
#' Cluster Name
#' 
#' @description Selects best fitting name for a cluster of pathways.
#' 
#' @param sim similarity matrix used to detect the clusters
#' @param clusters a list of clusters, where values are cluster ID and names are pathway names
#' @param method a method for setting cluster names. Available method include \code{'pagerank'},
#' \code{'hits'} and \code{'none'}
#' 
findClusterNames <- function(sim,
                             clusters,
                             method = c('pagerank', 'hits', 'none')) {
  method <- match.arg(method)

  switch(method,
         'pagerank' = clusterNamesPagerank(sim, clusters),
         'hits' = clusterNamesHits(sim, clusters),
         'none' = clusters
  )
}
