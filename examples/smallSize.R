# Testing with a small number of enrichment results

detach('package:pathExplore', unload = TRUE)
library(pathExplore)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)
data(geneList)

enrichment <- gseGO(geneList, 'CC', org.Hs.eg.db, eps = 0)
enrichment <- enrichment@result %>% setDT

# Select random pathways
set.seed(348934)
dt <- enrichment[ sample(1:nrow(enrichment), 10) ] %>% setDF

enrichmentNetwork(dt, drawEllipses = TRUE, clustMethod = 'spectral')
