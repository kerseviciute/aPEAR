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
dt <- enrichment[ sample(1:nrow(enrichment), 20) ] %>% setDF

# Example how to change colors
enrichmentNetwork(dt, drawEllipses = TRUE, clustMethod = 'hier') +
  scale_colour_gradient(low = '#0076A4', high = '#881D58')
