library(dplyr)
library(data.table)

message('Dataset 1')
dt <- readRDS('datasets/dataset1.RDS')
dt <- dt@result %>% setDT
dt <- dt[ order(pvalue) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/1.RDS')

message('Dataset 2')
dt <- readRDS('datasets/dataset2.RDS')
dt <- dt@result %>% setDT
dt <- dt[ order(pvalue) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/2.RDS')

message('Dataset 3')
dt <- readRDS('datasets/dataset3.RDS')
dt <- dt@result %>% setDT
dt <- dt[ order(pvalue) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/3.RDS')

message('Dataset 4')
dt <- readRDS('datasets/dataset4.RDS')
dt <- dt@result %>% setDT
dt <- dt[ order(pvalue) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/4.RDS')

message('Dataset 5')
dt <- readRDS('datasets/dataset5.RDS')
dt <- dt@result %>% setDT
dt <- dt[ order(pvalue) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/5.RDS')

message('Dataset 6')
dt <- readRDS('datasets/dataset6.RDS')
setnames(dt, c('Size', 'pathwayGenes'), c('setSize', 'core_enrichment'))
dt <- dt[ order(FDR) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/6.RDS')

message('Dataset 7')
dt <- readRDS('datasets/dataset7.RDS')
setnames(dt, c('Size', 'pathwayGenes'), c('setSize', 'core_enrichment'))
dt <- dt[ order(FDR) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/7.RDS')

message('Dataset 8')
dt <- readRDS('datasets/dataset8.RDS')
setnames(dt, c('Size', 'pathwayGenes'), c('setSize', 'core_enrichment'))
dt <- dt[ order(FDR) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/8.RDS')

message('Dataset 9')
dt <- readRDS('datasets/dataset9.RDS')
dt <- dt@result %>% setDT
dt <- dt[ order(pvalue) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/9.RDS')

message('Dataset 10')
dt <- readRDS('datasets/dataset10.RDS')
dt <- rbind(dt$DO@result, dt$KEGG@result) %>% setDT
dt <- dt[ order(pvalue) ][ 1:1000 ] %>% na.omit %>% setDF
stopifnot(nrow(dt) == 1000)
stopifnot(all(c('Description', 'core_enrichment', 'NES', 'setSize') %in% colnames(dt)))
saveRDS(dt, 'datasets/10.RDS')
