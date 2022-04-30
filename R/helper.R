emptyMatrix <- function(genes) {
  sim <- matrix(data = NA, nrow = length(genes), ncol = length(genes))
  colnames(sim) <- names(genes)
  rownames(sim) <- names(genes)

  sim
}

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