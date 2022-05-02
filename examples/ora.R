# Testing ORA analysis using clusterProfiler

# detach("package:pathExplore", unload=TRUE)
library(pathExplore)

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
data(geneList)

geneList <- geneList[ abs(geneList) > 2 ]

enrich <- enrichGO(names(geneList), OrgDb = org.Hs.eg.db, ont = 'MF')
enrich@result %>% nrow

p <- enrichmentNetwork(enrichment, verbose=T, colorType = 'pval', pCutoff = -5)
p
