
library(tidyverse)

in_f <- list.files('kaiju', full.names = TRUE)

data_in <- lapply(in_f, function(x){
  read_tsv(x) %>% mutate(id = sub(".out", "", file))
}) %>% bind_rows() %>% select(id, taxon_name, percent)

## mean 0.72 

cutoff = 0.01

heat <- data_in %>% group_by(taxon_name) %>% 
  filter(max(percent) >= cutoff) %>%
  mutate(detected = ifelse(percent >= cutoff, yes = 1, no = 0)) %>%
  select(-percent) %>%
  spread(taxon_name, detected) %>%
  column_to_rownames('id') %>%
  select(-contains('cannot'), -'unclassified')

heat

pheatmap::pheatmap(t(heat), cluster_rows = FALSE,
                   filename = 'output/kaiju_detection.pdf')
dev.off()
pheatmap::pheatmap(t(heat), cluster_rows = FALSE,
                   filename = 'output/kaiju_detection.png')
dev.off()


