library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

results <- readRDS('evaluation/simSpeed.RDS') %>%
	.[ , Method := factor(Method, levels = c('markov', 'hier', 'spectral', 'jaccard', 'cosine', 'cor'))] %>%
	.[ , Similarity := factor(Similarity, levels = c('jaccard', 'cosine', 'cor'))]

# time ~ similarity + N
# How long does it take to calculate the similarity matrix using method x and N data points?

# ggplot() +
# 	geom_line(data = dt, aes(x = N, y = Seconds, color = Method)) +
# 	geom_point(data = dt, aes(x = N, y = Seconds, color = Method)) +
# 	theme_bw()

# ggsave('similarity_time.png', plot)

# What is the time range at N data points for similarity method x?
plot1 <- results[ Type == 'similarity' ] %>%
	.[ , Seconds := Time ] %>%
	.[ , N := factor(Chunk) ] %>%
	.[ , Similarity := Method ] %>%
	ggplot(aes(x = N, y = Seconds, color = Similarity)) +
	geom_boxplot(position = 'identity', width = 0.5) +
	stat_summary(fun=mean, geom="line", aes(color = Similarity, group = Similarity)) +
	stat_summary(fun=mean, geom = 'point') +
	theme_bw() +
	theme(legend.position = 'bottom',
		legend.text=element_text(size=8), legend.title=element_text(size=9))

# ggsave('similarity_time_boxplot.png', plot)

# How long does it take to calculate clusters using methods x for N data points?
plot2 <- results[ Type == 'clustering'] %>%
	.[ , Seconds := Time ] %>%
	.[ , N := factor(Chunk) ] %>%
	.[ , Cluster := Method ] %>%
	ggplot(aes(x = N, y = Seconds, color = Cluster)) +
		geom_boxplot(position = 'identity', width = 0.5) +
		stat_summary(fun=mean, geom="line", aes(color = Cluster, group = Cluster)) +
		stat_summary(fun=mean, geom = 'point') +
		theme_bw() +
		theme(legend.position = 'bottom', legend.text=element_text(size=8), legend.title=element_text(size=9))

plot <- ggarrange(plot1, plot2, labels = c('A', 'B'), ncol = 2)

plot

ggsave('sim_clust_time.png', plot)

# ggsave('cluster_time.png', plot)

# Does time calculation depend on similarity method used?
plot3 <- results[ Type == 'clustering' ] %>%
	.[ , list(Seconds = mean(Time), N = Chunk), by = c('Chunk', 'Method', 'Similarity') ] %>%
	ggplot(aes(x = N, y = Seconds, color = Similarity)) +
	facet_wrap(Method ~ ., scales = 'free') +
	geom_line() +
	geom_point() +
	theme(legend.position = 'bottom') +
	theme_bw() +
	theme(legend.position = 'bottom', legend.text=element_text(size=8), legend.title=element_text(size=9))

plot3

ggsave('cluster_time_by_similarity.png', plot3)


results <- readRDS('evaluation/clusterScore.RDS') %>%
	.[ , Cluster := factor(Cluster, levels = c('markov', 'hier', 'spectral')) ] %>%
	.[ , Similarity := factor(Similarity, levels = c('jaccard', 'cosine', 'cor')) ]

# results %>%
#   reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
#                  variable.name = 'Score', value.name = 'Value') %>%
#   setDT %>%
#   ggplot(aes(x = Cluster, y = Value)) +
#   facet_wrap(Score ~ ., scales = 'free') +
#   geom_boxplot(outlier.alpha = 0) +
#   geom_jitter(width = 0.2, height = 0, aes(color = Similarity)) +
#   stat_compare_means(comparisons = list(c('markov', 'hier'),
#   	c('hier', 'spectral'), c('markov', 'spectral')), method = 'wilcox.test')

# results %>%
#   reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
#                  variable.name = 'Score', value.name = 'Value') %>%
#   setDT %>%
#   ggplot(aes(x = Similarity, y = Value)) +
#   facet_wrap(Score ~ ., scales = 'free') +
#   geom_boxplot(outlier.alpha = 0) +
#   geom_jitter(width = 0.2, height = 0, aes(color = Cluster)) +
#   stat_compare_means(comparisons = list(c('jaccard', 'cor'),
#   	c('jaccard', 'cosine'), c('cosine', 'cor')), method = 'wilcox.test')

plot1 <- results %>%
  reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
                 variable.name = 'Score', value.name = 'Value') %>%
  setDT %>%
  ggplot(aes(x = Cluster, y = Value)) +
  facet_wrap(Score ~ ., scales = 'free') +
  geom_boxplot(outlier.alpha = 0, width = 0.8, lwd = 0.2) +
  geom_boxplot(outlier.alpha = 0, aes(color = Similarity, fill = Similarity),
  	alpha = 0.2, lwd = 0.4,
  	position = position_dodge(width=0.6), width = 0.5) +
  stat_compare_means(comparisons = list(c('markov', 'hier'),
  	c('hier', 'spectral'), c('markov', 'spectral')), method = 'wilcox.test') +
  theme_bw() +
  theme(legend.position = 'bottom', axis.text.x = element_text(size = 8))

# ggsave('cluster_quality_compare.png', plot)

plot2 <- results %>%
  reshape2::melt(id.vars = c('Cluster', 'Similarity', 'Dataset'),
                 variable.name = 'Score', value.name = 'Value') %>%
  setDT %>%
  ggplot(aes(x = Similarity, y = Value)) +
  facet_wrap(Score ~ ., scales = 'free') +
  geom_boxplot(outlier.alpha = 0, width = 0.8, lwd = 0.2) +
  geom_boxplot(outlier.alpha = 0, aes(color = Cluster, fill = Cluster),
  	alpha = 0.2, lwd = 0.4,
	position = position_dodge(width=0.6), width = 0.5) +
  stat_compare_means(comparisons = list(c('jaccard', 'cosine'),
  	c('cosine', 'cor'), c('jaccard', 'cor')), method = 'wilcox.test') +
  theme_bw() +
  theme(legend.position = 'bottom', axis.text.x = element_text(size = 8))

plot <- ggarrange(plot1, plot2, labels = c('A', 'B'), ncol = 2)

plot

ggsave('quality_compare.png', plot)

results <- readRDS('evaluation/plotSpeed.RDS')

plot <- results %>%
    .[ , Seconds := mean(Time), by = c('Similarity', 'Cluster', 'Chunk') ] %>%
    .[ , N := Chunk ] %>%
    .[ , list(Similarity, Cluster, Seconds, N) ] %>%
    .[ , Similarity := factor(Similarity, levels = c('jaccard', 'cosine', 'cor')) ] %>%
    .[ , Cluster := factor(Cluster, levels = c('markov', 'hier', 'spectral')) ] %>%
    ggplot(aes(x = N, y = Seconds, shape = Cluster, color = Similarity,
        group=interaction(Similarity, Cluster))) +
    geom_line() +
    geom_point() +
    theme_bw()

plot

ggsave('plot_time.png', plot)
