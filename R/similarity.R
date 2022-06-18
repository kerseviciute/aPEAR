#' 
#' Pathway Jaccard Similarity
#' 
#' @description Calculates jaccard similarity between all pathways.
#' 
#' @param genes genes in all pathways. Should be a list of lists, which each entry containing the
#' genes that belong to a specific pathway (names of the list should be the pathways)
#' 
#' @importFrom bayesbio jaccardSets
#' 
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

#' 
#' Pathway Cosine Similarity
#' 
#' @description Calculates cosine similarity between all pathways.
#' 
#' @param genes genes in all pathways. Should be a list of lists, which each entry containing the
#' genes that belong to a specific pathway (names of the list should be the pathways)
#' 
#' @importFrom lsa cosine
#' 
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

#' 
#' Pathway Correlation
#' 
#' @description Calculates correlation between all pathways.
#' 
#' @param genes genes in all pathways. Should be a list of lists, which each entry containing the
#' genes that belong to a specific pathway (names of the list should be the pathways)
#' 
similarityCorrelation <- function(genes) {
  m <- occurenceMatrix(genes)
  sim <- cor(t(m))
  sim[ sim < 0 ] <- 0
  sim
}

#' 
#' Similarity Matrix
#' 
#' @description Calculates a similarity matrix of all pathways in the enrichment.
#' 
#' @param enrichment a data frame containing enrichment results
#' @param geneCol which column contains gene lists
#' @param pathCol which column contains path descriptions (human readable text, not ids!)
#' @param method a method to be used. Available values: \code{'jaccard'}, \code{'cosine'} and
#' \code{'cor'}
#' 
#' @importFrom tibble deframe
#' @importFrom dplyr %>%
#' 
#' @export
#' 
pathwaySimilarity <- function(
  enrichment,
  geneCol,
  pathCol = 'Description',
  method = c('jaccard', 'cosine', 'cor')
) {
  method <- match.arg(method)

  genes <- getGenes(enrichment, geneCol, pathCol)

  switch(method,
         'jaccard' = similarityJaccard(genes),
         'cosine' = similarityCosine(genes),
         'cor' = similarityCorrelation(genes)
  )
}

#'
#' Enrichment Genes
#'
#' @description Gets genes involved in each pathway of enrichment result.
#'
#' @param enrichment a data frame containing enrichment results
#' @param geneCol which column contains gene lists
#' @param pathCol which column contains path descriptions (human readable text, not ids!)
#'
#' @export
#'
getGenes <- function(enrichment, geneCol, pathCol = 'Description') {
  cols <- c(pathCol, geneCol)

  stopifnot(!is.null(geneCol))
  stopifnot(all(cols %in% colnames(enrichment)))

  enrichment[ , cols ] %>%
    deframe %>%
    lapply(function(x) strsplit(x, split = '/')[[1]])
}
