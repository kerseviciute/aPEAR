# TESTS ALL POSSIBLE COMBINATIONS

# detach("package:pathExplore", unload=TRUE)
library(pathExplore)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
data(geneList)
library(foreach)

enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'MF')

simMethods <- c('jaccard', 'cosine', 'cor')
clustMethods <- c('markov', 'hier', 'spectral')
namingMethods <- c('pagerank', 'hits')

enrichment <- enrich@result

for (simMethod in simMethods) {
  for (clustMethod in clustMethods) {
    for (clustNameMethod in namingMethods) {
      p <- enrichmentNetwork(enrich@result,
                             simMethod = simMethod,
                             clustMethod = clustMethod,
                             clustNameMethod = clustNameMethod,
                             verbose = FALSE
      )
      print(p)
    }
  }
}
