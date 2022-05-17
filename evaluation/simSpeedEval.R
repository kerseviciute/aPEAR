library(pathExplore)
library(foreach)
library(glue)
library(ggplot2)

setwd('evaluation')

chunks <- seq.int(100, 1000, 50)
datasets <- 1:10
methods <- c('jaccard', 'cosine', 'cor')

nCalculations <- length(chunks) * length(datasets) * length(methods)
message('Calculations awaiting: ', nCalculations)

set.seed(2395793)

results <- list()
counter <- 0
for (dataset in datasets) {
  dt <- readRDS(glue('datasets/{dataset}.RDS'))

  for (chunk in chunks) {
    for (method in methods) {
      counter <- counter + 1
      message(counter)

      message('Dataset ', dataset, ' chunk ', chunk, ' method ', method)
      time <- system.time(pathwaySimilarity(dt[ 1:chunk, ], geneCol = 'core_enrichment', method = method))
      time <- time[[ 'elapsed' ]]

      results[[ counter ]] <- data.table(Chunk = chunk, Dataset = dataset, Method = method, Time = time)
    }
  }
}

results <- rbindlist(results)

stopifnot(nrow(results) == nCalculations)

results[ , list(Seconds = mean(Time), N = Chunk), by = c('Chunk', 'Method') ] %>%
  ggplot(aes(x = N, y = Seconds, color = Method)) +
  geom_point() +
  geom_line()

saveRDS(results, 'simSpeed.RDS')
