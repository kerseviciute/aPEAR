library(pathExplore)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(foreach)
library(data.table)
data(geneList)

set.seed(45867)
enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'CC', eps = 0)
dim(enrich@result)

simMethods <- c('jaccard')
clustMethods <- c('markov')
namingMethods <- c('pagerank')

nCalculations <- length(simMethods) * length(clustMethods) * length(namingMethods)
message('Calculations awaiting: ', nCalculations)

counter <- 0
plots <- list()
data <- list()
for (simMethod in simMethods) {
  for (clustMethod in clustMethods) {
    for (clustNameMethod in namingMethods) {
      set.seed(349753)
      counter <- counter + 1
      plots[[ counter ]] <- enrichmentNetwork(enrich@result,
                             simMethod = simMethod,
                             clustMethod = clustMethod,
                             clustNameMethod = clustNameMethod,
                             drawEllipses = TRUE,
                             verbose = FALSE,
                             fontSize = 3,
                             innerCutoff = 0.005, outerCutoff = 0.2,
                             repelLabels = TRUE,
                             colorBy = 'p.adjust',
                             colorType = 'pval'
      )
      data[[ counter ]] <- data.table(Index = counter, Similarity = simMethod, Cluster = clustMethod, Name = clustNameMethod)

      message('Done ', counter, '/', nCalculations)
    }
  }
}

plots[[ counter ]]

data <- rbindlist(data)

# saveRDS(list(data = data, plots = plots), 'examplePlots.RDS')

for (i in seq_along(plots)) {
  filename <- paste0(paste(data[ i, Similarity ], data[ i, Cluster ], data[ i, Name ], sep = '_'), '.png')
  message('Saving to ', filename)
  ggsave(filename, plots[[ i ]])
}

enrichment <- readRDS('enrichment.RDS')
enrichment <- enrichment[ FDR < 0.02 ] %>% setDF

set.seed(285934)
plot <- enrichmentNetwork(enrichment,
                             simMethod = 'jaccard',
                             clustMethod = 'hier',
                             clustNameMethod = 'pagerank',
                             colorBy = 'NES',
                             nodeSize = 'Size',
                             drawEllipses = TRUE,
                             verbose = FALSE,
                             minClusterSize = 4,
                             fontSize = 3,
                             innerCutoff = 0.005, outerCutoff = 0.2,
                             repelLabels = TRUE)

ggsave('../images/enrichment_pathExplore.png', plot)




genes <- geneList[ abs(geneList) > 2 ] %>% names
enrich <- enrichGO(genes, OrgDb = org.Hs.eg.db, ont = 'BP')
dim(enrich@result)

library(enrichplot)
edo <- pairwise_termsim(enrich)
plot <- emapplot(edo, showCategory = 130, node_label = 'group', group_category = TRUE)
plot
ggsave('/Users/ieva/Thesis2022/pathExplore/images/emapplot.png', plot)

dt <- enrich@result
dt %>% setDT
dt <- dt[ qvalue < 0.01 ]
dt %>% setDF

plot <- enrichmentNetwork(dt,
                  simMethod = 'jaccard',
                  clustMethod = 'hier',
                  clustNameMethod = 'pagerank',
                  drawEllipses = TRUE,
                  verbose = TRUE,
                  fontSize = 3,
                  innerCutoff = 0.2, outerCutoff = 0.4,
                  repelLabels = TRUE, pCutoff = -15,
                  colorBy = 'p.adjust',
                  colorType = 'pval')
plot
ggsave('/Users/ieva/Thesis2022/pathExplore/images/pval_enrichmentNetwork.png', plot)







packages <- c(
'data.table',
'foreach',
'ggforce',
'ggplot2',
'ggrepel',
'igraph',
'MCL',
'arules',
'bayesbio',
'dplyr',
'lsa',
'tibble',
'plotly',
'clusterProfiler',
'org.Hs.eg.db',
'DOSE',
'glue',
'reshape2',
'clusterCrit',
'ggpubr'
)

library(glue)

for (package in packages) {
  print(glue('{package} v{packageVersion(package)}\\newline'))
}












