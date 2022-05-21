#
# Clustering evaluation.
# Evaluated for each clustering method using 500 pathways. Averaged for each similarity method
# and dataset.
#

pacman::p_load(data.table)
pacman::p_load(dplyr)
pacman::p_load(glue)
pacman::p_load(reshape2)
pacman::p_load(clusterCrit)
pacman::p_load(ggpubr)
library(pathExplore)

setwd('evaluation')

clustMethods <- list('markov', 'hier', 'spectral')
simMethods <- list('jaccard', 'cosine', 'cor')
datasets <- 1:10

nCalculations <- length(clustMethods) *
  length(simMethods) *
  length(datasets)
message('Calculations awaiting: ', nCalculations)

size <- 500

set.seed(248574)

result <- list()
counter <- 0
for (dataset in datasets) {
  dt <- readRDS(glue('datasets/{dataset}.RDS'))
  dt <- dt[ 1:size, ]


  for (simMethod in simMethods) {
    sim <- pathwaySimilarity(dt, geneCol = 'core_enrichment', method = simMethod)

    for (clustMethod in clustMethods) {
      message('Done ', counter, '/', nCalculations)
      message('Dataset ', dataset, ' method ', simMethod, ' clustMethod ', clustMethod)
      counter <- counter + 1

      clusters <- findClusters(sim, method = clustMethod, nameMethod = 'none')
      pathways <- names(clusters)
      clusters <- as.integer(clusters)

      idx <- colnames(sim) %in% pathways
      simClust <- sim[ idx, idx ]

      score1 <- intCriteria(simClust, clusters, 'Dunn')
      score2 <- intCriteria(simClust, clusters, 'Silhouette')
      score3 <- intCriteria(simClust, clusters, 'Davies_Bouldin')

      result[[ counter ]] <- data.table(Dataset = dataset,
                                        Cluster = clustMethod,
                                        Similarity = simMethod,
                                        Dunn = unlist(score1),
                                        Silhouette = unlist(score2),
                                        DaviesBouldin = unlist(score3))
    }
  }
}

results <- rbindlist(result)

results %>%
  reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
                 variable.name = 'Score', value.name = 'Value') %>%
  ggplot(aes(x = Cluster, y = Value)) +
  facet_grid(cols = vars(Score), space = 'free') +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, aes(color = Similarity))

results %>%
  reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
                 variable.name = 'Score', value.name = 'Value') %>%
  setDT %>%
  .[ Score == 'Dunn' ] %>%
  ggplot(aes(x = Similarity, y = Value)) +
  facet_grid(cols = vars(Cluster)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, aes(color = Cluster)) +
  stat_compare_means()

results %>%
  reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
                 variable.name = 'Score', value.name = 'Value') %>%
  setDT %>%
  .[ Score == 'Silhouette' ] %>%
  ggplot(aes(x = Cluster, y = Value)) +
  facet_grid(cols = vars(Similarity)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, aes(color = Similarity)) +
  stat_compare_means()

results %>%
  reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
                 variable.name = 'Score', value.name = 'Value') %>%
  setDT %>%
  .[ Score == 'DaviesBouldin' ] %>%
  ggplot(aes(x = Cluster, y = Value)) +
  facet_grid(cols = vars(Similarity)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, aes(color = Similarity)) +
  stat_compare_means()

saveRDS(results, 'clusterScore.RDS')
