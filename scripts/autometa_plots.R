# Load packages
library("ggplot2")

# Read the main binning table

f_in <- list.files('Autometa', full.names = TRUE)

#filepath="Autometa/AttbisABBM3_FD.bacteria.hdbscan.main.tsv"


data_l = lapply(f_in, function(x){
  d_in <- read.table(x, header=TRUE, sep='\t')
  d_in$id <- sub(".bacteria.hdbscan.main.tsv", "", basename(x))
  d_in
})

data <- do.call(rbind, data_l)


data_filtered <- data %>% group_by(order) %>%
  mutate(occurence = length(unique(id))) %>%
  ungroup() %>%
  filter(occurence >= 6)

# Fill empty cells as unclustered
data_filtered$cluster <- sub("^$", "Unclustered", data$cluster)

ggplot(data_filtered, aes(x=x_1, y=x_2, color=order, group=order)) +
  geom_point(size=(sqrt(data_filtered$length))/100, shape=20, alpha=0.5) +
  theme_classic() + 
  xlab('BH-tSNE X') + 
  ylab('BH-tSNE Y') +
  guides(color = guide_legend( title = 'Bacterial Order')) +
  facet_wrap(~id) +
  ggsave(filename = 'output/autometa.png')

ggplot(data_filtered, aes(x=x_1, y=x_2, color=cluster, group=cluster)) +
  geom_point(size=(sqrt(data_filtered$length))/100, shape=20, alpha=0.5) +
  theme_classic() + 
  xlab('BH-tSNE X') + 
  ylab('BH-tSNE Y') +
  guides(color = guide_legend( title = 'Bacterial Clusters')) +
  facet_wrap(~id) +
  ggsave(filename = 'output/autometa_bins.png')
          