library(pathExplore)
library(foreach)
library(glue)
library(ggplot2)

setwd('evaluation')

chunks <- seq.int(100, 500, 50)
datasets <- 1:10
simMethods <- c('jaccard', 'cosine', 'cor')
clustMethods <- c('markov', 'hier', 'spectral')

nCalculations <- length(chunks) * length(datasets) * length(simMethods) +
  length(chunks) * length(datasets) * length(simMethods) * length(clustMethods)
message('Calculations awaiting: ', nCalculations)

set.seed(2395793)

results <- list()
counter <- 0
for (dataset in datasets) {
  dt <- readRDS(glue('datasets/{dataset}.RDS'))

  for (chunk in chunks) {
    for (simMethod in simMethods) {
      message('Done ', counter, '/', nCalculations)
      counter <- counter + 1

      message('Dataset ', dataset, ' chunk ', chunk, ' method ', simMethod)
      time <- system.time({
        sim <- pathwaySimilarity(dt[ 1:chunk, ], geneCol = 'core_enrichment', method = simMethod)
      })
      simTime <- time[[ 'user.self' ]]

      results[[ counter ]] <- data.table(Chunk = chunk,
                                         Dataset = dataset,
                                         Method = simMethod,
                                         Time = simTime,
                                         Type = 'similarity',
                                         Similarity = NA
      )

      for (clustMethod in clustMethods) {
        set.seed(2395029)

        counter <- counter + 1
        message('Dataset ', dataset, ' chunk ', chunk, ' method ', simMethod, ' clustMethod ', clustMethod)

        # Sometimes spectral clustering fails due to randomization in kmeans algorithm, so
        # we repeat until it doesnt fail
        clusters <- NULL
        while (is.null(clusters)) {
          time <- system.time({
            clusters <- tryCatch({
              findClusters(sim, method = clustMethod, nameMethod = 'none')
            }, error = function(x) {
              NULL
            })
          })
        }

        clustTime <- time[[ 'user.self' ]]

        results[[ counter ]] <- data.table(Chunk = chunk,
                                           Dataset = dataset,
                                           Method = clustMethod,
                                           Time = clustTime,
                                           Type = 'clustering',
                                           Similarity = simMethod
        )
      }
    }
  }
}

results <- rbindlist(results)

stopifnot(nrow(results) == nCalculations)

results[ Type == 'clustering'] %>%
  .[ Similarity == 'cor' ] %>%
  .[, list(Seconds = mean(Time), N = Chunk), by = c('Chunk', 'Method') ] %>%
  ggplot(aes(x = N, y = Seconds, color = Method)) +
  geom_point() +
  geom_line()


results[ Type == 'similarity', list(Seconds = mean(Time), N = Chunk), by = c('Chunk', 'Method') ] %>%
  ggplot(aes(x = N, y = Seconds, color = Method)) +
  geom_point() +
  geom_line()

saveRDS(results, 'simSpeed.RDS')
