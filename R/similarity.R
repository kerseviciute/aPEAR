#' Calculates jaccard similarity between all pathways.
#'
#' @importFrom bayesbio jaccardSets
similarityJaccard <- function(genes) {
  sim <- emptyMatrix(genes)

  for (i in seq_along(genes)) {
    for (j in seq_along(genes)) {
      if (j > i) { break }
      sim[ i, j ] <- jaccardSets(genes[[i]], genes[[j]])
      sim[ j, i ] <- sim[ i, j ]
    }
  }

  sim
}

#' Calculates cosine similarity between all pathways.
#'
#' @importFrom lsa cosine
similarityCosine <- function(genes) {
  sim <- emptyMatrix(genes)
  m <- occurenceMatrix(genes)

  for (i in seq_along(genes)) {
    sim[ i, i ] <- 1.0
    for (j in seq_along(genes)) {
      if (j >= i) { break }
      terms <- c(genes[[i]], genes[[j]]) %>% unique
      t <- m[ c(i, j), colnames(m) %in% terms ]
      sim[ i, j ] <- cosine(t[ 1, ], t[ 2, ])
      sim[ j, i ] <- sim[ i, j ]
    }
  }

  sim
}

#' Calculates mutual information between all pathways.
#'
#' @importFrom infotheo mutinformation
similarityMutualInformation <- function(genes) {
  sim <- emptyMatrix(genes)
  m <- occurenceMatrix(genes)

  for (i in seq_along(genes)) {
    for (j in seq_along(genes)) {
      if (j > i) { break }
      sim[ i, j ] <- mutinformation(m[ i, ], m[ j, ], method = 'emp')
      sim[ j, i ] <- sim[ i, j ]
    }
  }

  sim
}

#' Calculates a similarity matrix.
#'
#' @description Calculates a similarity matrix of all pathways in the enrichment.
#'
#' @param enrichment a data frame containing enrichment results.
#' @param method a method to be used. Available values: 'jaccard', 'cosine', 'mutualInfo'.
#'
#' @importFrom tibble deframe
#' @importFrom dplyr %>%
#'
#' @export
pathwaySimilarity <- function(enrichment, method = 'jaccard') {
  availableSimilarityMethods <- c('jaccard', 'cosine', 'mutualInfo')
  if (!(method %in% availableSimilarityMethods)) {
    stop(paste0('Unavailable method "', method, '"'))
  }

  # TODO: how to add the option to use all genes in the pathway?
  # cols <- c('Description', ifelse(useOnlyCore, 'core_enrichment', 'enrichment'))
  cols <- c('Description', 'core_enrichment')
  genes <- enrichment[ , cols ] %>%
    deframe %>%
    lapply(\(x) strsplit(x, split = '/')[[1]])

  if (method == 'jaccard') return(similarityJaccard(genes))
  if (method == 'cosine') return(similarityCosine(genes))
  if (method == 'mutualInfo') return(similarityMutualInformation(genes))
}
