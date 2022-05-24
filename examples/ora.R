# Testing ORA analysis using clusterProfiler

# detach("package:pathExplore", unload=TRUE)
library(pathExplore)

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(data.table)
library(dplyr)
data(geneList)

geneList <- geneList[ abs(geneList) > 2 ]

enrich <- enrichGO(names(geneList), OrgDb = org.Hs.eg.db, ont = 'MF')
enrichment <- enrich@result %>%
  setDT %>%
  .[ order(pvalue) ] %>%
  .[ 1:200 ] %>%
  setDF

p <- enrichmentNetwork(enrichment,
                       clustMethod = 'hier',
                       verbose = T,
                       colorType = 'pval',
                       pCutoff = -5,
                       fontSize = 2,
                       drawEllipses = TRUE,
                       minClusterSize = 4
)
p
