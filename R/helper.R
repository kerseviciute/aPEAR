#'
#' Creating Valid Input
#' 
#' @description Currently, this package is able to work with DOSE enrichment classes. However, you
#' may wish to pass different enrichment results. In this case, prepare your data so that it
#' contains \code{'Description'} (readable pathway names) and \code{'pathwayGenes'} (a list of genes
#' that belong to the pathway, separated by '/') columns. Use \code{'colorBy'} and \code{'nodeSize'}
#' parameters to set which columns should be used for colouring the nodes and adjusting their size
#' (eg. \code{colorBy = 'NES', nodeSize = 'setSize'}). These will not be set automatically if you
#' create custom input. You may use \code{validateEnrichment} to check whether your data has valid
#' format.
#' 
#' @seealso \code{enrichmentNetwork}, \code{validateEnrichment}
#' @name enrichmentData
#' 
NULL

#'
#' Empty Matrix
#' 
#' @description Creates an empty matrix from names of a list. Will set colnames and rownames as
#' the names of the list.
#' 
#' @param values a list of values
#' @param data default data
#' 
emptyMatrix <- function(values, data = NA) {
  sim <- matrix(data = data, nrow = length(values), ncol = length(values))
  colnames(sim) <- make.unique(names(values))
  rownames(sim) <- make.unique(names(values))

  sim
}

#' 
#' Occurence Matrix
#' 
#' @description Creates an occurence matrix from gene lists (rows - pathways, cols - genes).
#' 
#' @param values a list of lists
#' 
#' @importFrom dplyr %>%
#' 
occurenceMatrix <- function(values) {
  uniqueGenes <- values %>% unlist %>% unique

  m <- matrix(data = FALSE, nrow = length(values), ncol = length(uniqueGenes))
  rownames(m) <- names(values)
  colnames(m) <- uniqueGenes

  for (i in seq_along(values)) {
    m[ i, colnames(m) %in% values[[ i ]] ] <- TRUE
  }

  m
}


#' 
#' Enrichment Validation
#' 
#' @description Validates enrichment input for use with \code{enrichmentNetwork} method and checks
#' for available columns. Currently only validates clusterProfiler results.
#' 
#' @param enrichment data frame containing enrichment results
#' @param colorBy which column will be used to color \code{enrichmentNetwork} nodes. Will try to set
#' automatically if it is NULL
#' @param nodeSize which column will be used to set size of \code{enrichmentNetwork} nodes. Will try
#' to set automatically if it is NULL
#' @param verbose enable / disable log messages
#'
#' @seealso \code{enrichmentData}
#' 
#' @export
#' 
validateEnrichment <- function(enrichment,
                               colorBy = NULL,
                               nodeSize = NULL,
                               verbose = TRUE
) {
  enrichmentType <- NULL

  if (all(c('Description', 'pathwayGenes') %in% colnames(enrichment))) {
    enrichmentType <- 'custom'
    genesCol <- 'pathwayGenes'

    stopifnot(!any(c(is.null(colorBy), is.null(nodeSize))))

    if (verbose) message('Setting colorBy to ', colorBy)
    if (verbose) message('Setting colorBy to ', colorBy)
  }

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
  }

  if (!(colorBy %in% colnames(enrichment))) {
    stop('"', colorBy, '" not found in enrichment columns.')
  }

  # Setting nodeSize
  if (is.null(nodeSize)) {
    if (enrichmentType == 'enrichDOSE') nodeSize <- 'Count'
    if (enrichmentType == 'gseDOSE') nodeSize <- 'setSize'

    if (verbose) message('Setting nodeSize to ', nodeSize, ' automatically')
  }

  if (!(nodeSize %in% colnames(enrichment))) {
    stop('"', nodeSize, '" not found in enrichment columns.')
  }

  return(list(
    class = enrichmentType,
    colorBy = colorBy,
    nodeSize = nodeSize,
    genesCol = genesCol
  ))
}

#'
#' Word Splitting
#'
#' @description Splits pathway name into multiple lines.
#'
#' @param x string to split
#'
splitWords <- function(x) {
  strwrap(x, width = 30) %>% paste(collapse = '\n')
}
