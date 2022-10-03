library(tidyverse)


files_in <- list.files('data', recursive = TRUE, pattern = 'card', 
                       full.names = TRUE)



data_in <- lapply(files_in, function(x){
  read_tsv(x) %>% 
    select(ORF_ID, Cut_Off, Best_Hit_ARO, `AMR Gene Family`) %>%
    mutate(sample = basename(dirname(x)))
}) %>% bind_rows()
