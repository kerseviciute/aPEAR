# pathExplore - an R package for high level GSEA exploration

![Alt text](results/example.png "Pathway cluster network")

## Installation

```R
library(devtools)
install_github('ievaKer/pathExplore')
```

### Run an example

```R
library(pathExplore)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
data(geneList)

enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'CC')
p <- enrichmentNetwork(enrich@result)

p
```

### Plotly integration

To create interactive plots, use plotly:

```R
library(plotly)

ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
```

You may find example plots in __results__ directory of this repository.
