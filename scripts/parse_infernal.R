library(tidyverse)


files_in <- list.files('data/', recursive = TRUE, pattern = 'infernal', 
                       full.names = TRUE)



data_in <- lapply(files_in, function(x){
  df <- read.delim(x, sep = "" , header = F ,
             na.strings ="", stringsAsFactors= F, skip = 2) %>% 
    select('V2', 'V3', 'V4', 'V6', 'V8', 'V9', 'V10', 'V11', 'V17', 'V18') %>%
    na.omit()
  colnames(df) = c('desc', 'rna_fam', 'target_scaffold', 
                   'clan',
                   'model_start',
                   'model_end',
                   'seq_start', 'seq_end',
                   'score', 'evalue')
  
  df %>%
    group_by(clan, rna_fam) %>%
    mutate(locus = paste(target_scaffold, seq_start, seq_end)) %>%
    summarize(num_hits = length(unique(locus))) %>%
    ungroup() %>%
    unite(clan, rna_fam, col = 'id', sep = ":") %>%
    mutate(id = sub("-:", "", id)) %>%
    mutate(sample = sub("_infernal-genome.tblout", "", basename(x)))
  
}) %>% bind_rows()


df_wide <- data_in %>% spread(sample, num_hits) %>%
  column_to_rownames('id')

df_wide[is.na(df_wide)] <- 0

df_wide_filt <- df_wide[which(rowSums(df_wide) > 10),]


df_wide_filt[df_wide_filt >= 10] <- 10


library(pheatmap)

pheatmap(df_wide_filt)






