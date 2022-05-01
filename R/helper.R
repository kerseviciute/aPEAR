#' Creating Valid Input
#'
#' @description Currently, this package is able to work with DOSE enrichment classes. However, you
#' may wish to pass different enrichment results. In this case, prepare your data so that it
#' contains 'Description' and 'core_enrichment' columns. Use \code{'colorBy'} and \code{'nodeSize'}
#' parameters to set which columns should be used for colouring the nodes and adjusting their size
#' (eg. \code{colorBy = 'NES', nodeSize = 'setSize'}).
#'
#' @seealso \code{enrichmentNetwork}
#' @name enrichmentData
NULL

#' Creates an empty matrix (with NA vals) from gene lists.
emptyMatrix <- function(genes, data = NA) {
  sim <- matrix(data = data, nrow = length(genes), ncol = length(genes))
  colnames(sim) <- names(genes)
  rownames(sim) <- names(genes)

  sim
}

#' Creates an occurence matrix from gene lists.
#'
#' @importFrom dplyr %>%
occurenceMatrix <- function(genes) {
  uniqueGenes <- genes %>% unlist %>% unique

  m <- matrix(data = FALSE, nrow = length(genes), ncol = length(uniqueGenes))
  rownames(m) <- names(genes)
  colnames(m) <- uniqueGenes

  for (path in rownames(m)) {
    m[ path, colnames(m) %in% genes[[path]] ] <- TRUE
  }

  m
}

validateEnrichment <- function(enrichment, colorBy = NULL, nodeSize = NULL, verbose = TRUE) {
  enrichmentType <- NULL

  if (all(c('Description', 'geneID', 'pvalue', 'Count') %in% colnames(enrichment))) {
    enrichmentType <- 'enrichDOSE'
    genesCol <- 'geneID'
  }

  if (all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(enrichment))) {
    enrichmentType <- 'gseDOSE'
    genesCol <- 'core_enrichment'
  }

  if (is.null(enrichmentType)) {
    stop('Unrecognized enrichment type.')
  }

  if (verbose) message('Detected enrichment type ', enrichmentType)

  # Setting colorBy
  if (is.null(colorBy)) {
    if (enrichmentType == 'enrichDOSE') colorBy <- 'pvalue'
    if (enrichmentType == 'gseDOSE') colorBy <- 'NES'

    if (verbose) message('Setting colorBy to ', colorBy, ' automatically')
  } else {
    if (!(colorBy %in% colnames(enrichment))) {
      stop('"', colorBy, '" not found in enrichment columns.')
    }
  }

  # Setting nodeSize
  if (is.null(nodeSize)) {
    if (enrichmentType == 'enrichDOSE') nodeSize <- 'Count'
    if (enrichmentType == 'gseDOSE') nodeSize <- 'setSize'

    if (verbose) message('Setting nodeSize to ', nodeSize, ' automatically')
  } else {
    if (!(nodeSize %in% colnames(enrichment))) {
      stop('"', nodeSize, '" not found in enrichment columns.')
    }
  }

  return(list(
    class = enrichmentType,
    colorBy = colorBy,
    nodeSize = nodeSize,
    genesCol = genesCol
  ))
}
