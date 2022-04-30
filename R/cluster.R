getClusters <- function(mSim, minClusterSize = 2, method = 'wordcloud') {
  if (minClusterSize < 1) {
    stop(paste0('minClusterSize = ', minClusterSize, '. Please choose minClusterSize >= 1.'))
  }

  pacman::p_load(MCL)

  allow1 <- ifelse(minClusterSize == 1, TRUE, FALSE)

  res <- mcl(mSim,
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

getClustersHierarchical <- function(sim, minClusterSize=2, method = 'wordcloud') {
  hCluster <- as.dist(1 - sim) %>% hclust
  clusters <- cutree(hCluster, h=0.9)

  if (minClusterSize != 1) {
    sizes <- clusters %>% table
    ids <- names(sizes)[ sizes >= minClusterSize ]
    clusters <- clusters[ clusters %in% ids ]
  }

  setClusterNames(sim, clusters, method)
}

getClustersSpectral <- function(sim, method = 'wordcloud') {
  pacman::p_load(Spectrum)
  clusters <- Spectrum(sim, maxk=50, showres=FALSE)
  clusters <- clusters$assignments
  names(clusters) <- colnames(sim)

  setClusterNames(sim, clusters, method)
}

setClusterNames <- function(data, clusters, method = 'pagerank') {
  methods <- c('none', 'pagerank', 'wordcloud', 'hits')
  if (!(method %in% methods)) {
    stop(paste0('Unavailable method "', method, '"'))
  }

  library(foreach)

  if (method == 'pagerank')  { return(clusterNamesPagerank(data, clusters)) }
  if (method == 'wordcloud') { return(clusterNamesWordcloud(clusters)) }
  if (method == 'hits')      { return(clusterNamesHits(data, clusters))}
  if (method == 'none')      { return(clusters) }
}

clusterNamesWordcloud <- function(clusters) {
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

clusterNamesPagerank <- function(data, clusters) {
  pacman::p_load(igraph)

  stopifnot(rownames(data) == colnames(data))

  paths <- rownames(data)
  edges <- list()
  counter <- 1

  for (i in 1:(nrow(data) - 1)) {
    for (j in (i + 1):ncol(data)) {
      value <- data[ i, j ]

      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        edges[[counter]] <- data.table(from=paths[ i ], to=paths[ j ], weight=value)
        counter <- counter + 1
      }
    }
  }

  edges <- rbindlist(edges)
  g <- graph_from_data_frame(edges, directed = FALSE)
  scores <- page.rank(g)$vector

  clusterNames <- foreach(cluster = unique(clusters), .combine = rbind) %do% {
    name <- scores[ names(scores) %in% names(clusters[ clusters == cluster ]) ] %>%
        which.max %>% names

    data.table(ClusterID = cluster, Name = name)
  }

  as.data.table(clusters, keep.rownames = TRUE) %>%
      merge(clusterNames, by.x = 'clusters', by.y = 'ClusterID') %>%
      .[ , list(rn, Name) ] %>%
      deframe
}

clusterNamesHits <- function(data, clusters) {
  pacman::p_load(arules)

  stopifnot(rownames(data) == colnames(data))
  adjacency <- matrix(data = 0, nrow = nrow(data), ncol = ncol(data))
  rownames(adjacency) <- rownames(data)
  colnames(adjacency) <- colnames(data)

  paths <- rownames(data)

  for (i in 1:(nrow(data) - 1)) {
    for (j in (i + 1):ncol(data)) {
      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        adjacency[ i, j ] <- 1
        adjacency[ j, i ] <- 1
      }
    }
  }

  scores <- arules::hits(adjacency, type = 'relative')

  clusterNames <- foreach(cluster = unique(clusters), .combine = rbind) %do% {
    name <- scores[ names(scores) %in% names(clusters[ clusters == cluster ]) ] %>%
        which.max %>% names

    data.table(ClusterID = cluster, Name = name)
  }

  as.data.table(clusters, keep.rownames = TRUE) %>%
      merge(clusterNames, by.x = 'clusters', by.y = 'ClusterID') %>%
      .[ , list(rn, Name) ] %>%
      deframe
}
