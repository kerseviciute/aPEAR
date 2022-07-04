# Testing with a small number of enrichment results

detach('package:pathExplore', unload = TRUE)
library(pathExplore)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)
library(plotly)
data(geneList)

enrichment <- gseGO(geneList, 'CC', org.Hs.eg.db, eps = 0)
enrichment <- enrichment@result %>% setDT

# Select random pathways
set.seed(348934)
dt <- enrichment %>% setDF

# Example how to change colors
enrichmentNetwork(dt,
                  drawEllipses = TRUE,
                  clustMethod = 'hier',
                  repelLabels = FALSE,
                  fontSize = 3,
                  outerCutoff = 0.3) +
  scale_colour_gradient(limits = c(-2.8, 2.8), low = '#0076A4', high = '#881D58')
