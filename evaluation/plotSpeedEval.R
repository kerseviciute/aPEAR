library(pathExplore)
library(foreach)
library(glue)
library(ggplot2)
library(data.table)

setwd('evaluation')

chunks <- seq.int(100, 500, 50)
datasets <- 1:10
simMethods <- c('jaccard', 'cosine', 'cor')
clustMethods <- c('markov', 'hier', 'spectral')

nCalculations <- length(chunks) *
    length(datasets) *
    length(simMethods) *
    length(clustMethods)
message('Calculations awaiting: ', nCalculations)

set.seed(2395793)

results <- list()
counter <- 0
for (dataset in datasets) {
  dt <- readRDS(glue('datasets/{dataset}.RDS'))

  for (chunk in chunks) {
    for (simMethod in simMethods) {
      for (clustMethod in clustMethods) {
        set.seed(2395029)
        message('Dataset ', dataset, ' chunk ', chunk, ' method ', simMethod, ' clustMethod ', clustMethod)

        # Sometimes spectral clustering fails due to randomization in kmeans algorithm, so
        # we repeat until it doesnt fail
        plot <- NULL
        while (is.null(plot)) {
          time <- system.time({
            plot <- tryCatch({
              enrichmentNetwork(dt[ 1:chunk, ],
                                simMethod = simMethod,
                                clustMethod = clustMethod,
                                clustNameMethod = 'none',
                                verbose = FALSE,
                                drawEllipses = FALSE
              )
            }, error = function(x) {
              message(x)
              message('Error occured, repeating...')
              NULL
            })
          })
        }

        time <- time[[ 'user.self' ]]

        counter <- counter + 1
        results[[ counter ]] <- data.table(Chunk = chunk,
                                           Dataset = dataset,
                                           Cluster = clustMethod,
                                           Similarity = simMethod,
                                           Time = time
        )

        message('Done ', counter, '/', nCalculations)
      }
    }
  }
}

results <- rbindlist(results)

stopifnot(nrow(results) == nCalculations)

saveRDS(results, 'plotSpeed.RDS')
