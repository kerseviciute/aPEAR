# _aPEAR_: an R package for enrichment network visualisation

![Alt text](example.png "Pathway cluster network")

## Installation

### Install latest stable version

```R
library(devtools)
install_github('ievaKer/aPEAR')
```

### Install latest beta version

```R
library(devtools)
install_github('ievaKer/aPEAR@development')
```

### Run an example

```R
library(aPEAR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
data(geneList)

enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'CC')
p <- enrichmentNetwork(enrich@result, drawEllipses = TRUE, fontSize = 2.5)

p
```

### Plotly integration

To create interactive plots, use plotly:

```R
library(plotly)

ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
```

You may find example plots in __results__ directory of this repository.
