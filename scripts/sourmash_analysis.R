library(tidyverse)
library(pheatmap) 

files_read <- list.files('sourmash/taxa_reads', pattern = 'csv', full.names = TRUE)
files_ass <- list.files('sourmash/genome_genome', pattern = 'csv', full.names = TRUE)



# selecting specified query k=51
# loaded query: 11231.5.197616.GTCCGC.fastq.gz... (k=51, DNA)
# loaded 82 signatures.
# 
# Starting prefetch sweep across databases.
# Found 42 signatures via prefetch; now doing gather.
# 
# overlap     p_query p_match
# ---------   ------- -------
#   3.9 Mbp        0.3%   91.1%    GCA_002928325.1.fna
# 3.6 Mbp        0.2%   73.5%    GCF_011044475.1.fna
# 2.9 Mbp        0.2%   54.5%    GCA_002928355.1.fna
# 1.5 Mbp        0.1%   42.7%    GCA_011753235.1.fna
# 2.8 Mbp        0.1%   17.4%    GCA_011752635.1.fna
# 2.9 Mbp        0.0%    8.5%    GCA_002928365.1.fna
# 3.6 Mbp        0.0%   11.1%    GCA_011752675.1.fna
# 430.0 kbp      0.0%    8.7%    GCA_011752745.1.fna
# 340.0 kbp      0.0%    0.3%    GCA_022457215.1.fna
# 300.0 kbp      0.0%    4.3%    GCA_011752845.1.fna
# 2.9 Mbp        0.0%    3.1%    GCA_011752515.1.fna
# 130.0 kbp      0.0%    1.8%    GCA_011752995.1.fna
# 50.0 kbp       0.0%    0.6%    GCA_011752985.1.fna
# 80.0 kbp       0.0%    0.4%    GCA_011754595.1.fna
# 60.0 kbp       0.0%    0.2%    GCA_011752625.1.fna
# 70.0 kbp       0.0%    0.2%    GCA_011752685.1.fna
# 60.0 kbp       0.0%    0.2%    GCA_011752715.1.fna
# 40.0 kbp       0.0%    0.2%    GCA_011752835.1.fna
# 20.0 kbp       0.0%    0.1%    GCA_011753085.1.fna
# 10.0 kbp       0.0%    0.2%    GCA_011753195.1.fna
# 10.0 kbp       0.0%    0.0%    ICBG2046.fna
# 10.0 kbp       0.0%    0.0%    ICBG736.fna
# found less than 1 bp  in common. => exiting
# 
# found 22 matches total;
# the recovered matches hit 1.0% of the query (unweighted)

## intersect_bp (overlap between genome and reads)
## f_orig_query (p_query) (how much of the original query belongs to the match. )
## f_match (p_match) how much of the match is in the remaining query metagenome,
## f_match_orig  is how much of the match is in the original query.
## f_unique_to_query, is the fraction of the database match that is unique with respect to the original query. 


t = 'acc a long sent'


join_type_acc = function(x){
  parts = strsplit(x, ' ')[[1]]
  paste0(parts[1], ' ', paste0(parts[-1], collapse = '-'))
}



data_reads<- lapply(files_read, function(x){
  id = paste0(sub("\\.sbt", "", basename(x)))
  read = strsplit(id, "_")[[1]][1]
  db = strsplit(id, "_")[[1]][2]
  read_csv(x, col_types = 'nnnnnnnncccnnnncccn') %>%
    mutate(id = read ) %>%
    mutate(db = db) %>%
    group_by(name) %>%
    mutate(name = ifelse(grepl('Type', name), 
                         yes = join_type_acc(name),
                         no = name)) %>%
    mutate(name = ifelse(grepl('phage', name), 
                         yes = join_type_acc(name),
                         no = name)) %>%
    mutate(name = ifelse(grepl("Bathyarchaeota", name), 
                         yes = join_type_acc(name),
                         no = name)) %>%
    mutate(name = ifelse(grepl("\\[", name), 
                         yes = sub("\\]", "", sub("\\[", "", name)),
                         no = name)) %>%
    ungroup()
} ) %>% bind_rows() %>%
  select(id, name, db, query_filename, intersect_bp, f_orig_query, f_match, query_bp) %>%
  separate(name, into = c('acc', 'genus', 'species', 'strain', 'etc'), 
            sep = ' ', remove = FALSE) %>%
  mutate(genus = ifelse(genus == 'uncultured', yes = paste0(genus, "_", species),
                        no = genus))


write_csv(data_reads, 'output/sourmash_reads_taxa.csv')

# data_ass <- lapply(files_ass, function(x){
#   read_csv(x, col_types = 'nnnnnnnncccnnnncccn') %>%
#     mutate(id = paste0(sub("\\..*$", "", basename(x)), "_assembly"))
# } ) %>% bind_rows() %>%
#   select(id, name, query_filename, intersect_bp, f_orig_query, f_match, query_bp)
# 
data_in <- bind_rows(data_reads)#, data_ass)
# 
# id_map <- data_in %>% select(id, query_filename) %>% unique()
# 
# data_in <- data_in %>% select(-query_filename)
# 
# thres_1 <- data_in %>%
#   ggplot(aes( x = f_match)) + geom_histogram(binwidth = .01) +
#   theme_classic() +
#   geom_vline(xintercept = .015)
# 
# thres_2 <- data_in %>%
#   ggplot(aes( x = intersect_bp)) + geom_histogram(bins = 100) +
#   theme_classic() +
#   geom_vline(xintercept = 30000)
# 
# library(cowplot)
# 
# cowplot::plot_grid(thres_1, thres_2) %>%
#   cowplot::save_plot(filename = 'output/cutoffs.pdf')
# 
# genomes_surv <- read_csv('tables/genomes_surveyed.csv')


#ids <- paste(genomes_surv$organism, genomes_surv$genbank, sep = ":")
#names(ids) <- genomes_surv$genbank

filt_data <- data_in %>% 
#  filter(f_match >= .015, intersect_bp >= 30000) %>%
  mutate(detected = 1) %>%
  mutate(name = sub('.fna', '', name)) %>%
  select(id, genus, db, detected) %>%
  group_by(db, id, genus) %>%
  summarize(detected = 1 ) %>%
  ungroup() %>%
  group_by(db, genus) %>%
  filter(sum(detected) > 2) %>%
  ungroup() %>% 
  spread(id, detected) %>% 
  filter(db != 'genbank-2022.08-plant-k31.csv')




db_map <- select(filt_data, db, genus) %>% distinct() %>%
  column_to_rownames('genus')




#filt_data$row_id <- ids[filt_data$name]

# filt_data <- filt_data %>% separate(row_id, into = c('organism', 'gbk'), sep = ':') %>%
#   mutate(genus = sub(" .*$", "", organism)) %>%
#   group_by(genus) %>%
#   mutate(detected = ifelse(1 %in% detected, yes = 1, no = 0)) %>%
#   ungroup() %>%
#   select(id, genus, detected) %>%
#   unique() %>%
#   spread(id, detected) %>%
#   column_to_rownames('genus')

filt_data <- filt_data %>% select(-db) %>%
  column_to_rownames('genus')


filt_data[is.na(filt_data)] = 0

# row_order <-  c("Leucoagaricus gongylophorus Ac12 GCA_000382605.1",
#     "Leucoagaricus gongylophorus GCA_022457215.1",
#     "Clavispora lusitaniae GCF_000003835.1_ASM383v1_genomic",
#     "Pantoea GCA_002928325.1",
#     "Pantoea GCA_002928355.1",
#     "Pantoea GCA_002928365.1",
#     "Pantoea GCA_011752515.1",
#     "Pantoea GCA_011752615.1",
#     "Pantoea GCA_011752625.1",
#     "Pantoea GCA_011752635.1",
#     "Pantoea GCA_011752675.1",
#     "Pantoea GCA_011752685.1",
#     "Pantoea GCA_011753045.1",
#     "Pantoea stewartii GCF_011044475.1",
#     "Enterobacter GCA_011752745.1",
#     "Enterobacter GCA_011752765.1",
#     "Klebsiella GCA_011752715.1",
#     "Klebsiella GCA_011752755.1",
#     "Cedecea GCA_011752845.1",
#     "Burkholderia GCA_011752995.1",
#     "Burkholderia GCA_011753085.1",
#     "Burkholderia GCA_011753175.1",
#     "Asaia GCA_011753235.1",
#     "Acinetobacter GCA_011753255.1",
#     "Klebsiella variicola F2R9 GCF_009648975.1"
#   )


# cols_ord <- c(
#   "AttbisABBM1_FD_assembly",
#   "AttbisABBM1_FD_reads",
#   "AttbisABBM2_FD_assembly",
#   "AttbisABBM2_FD_reads",
#   "AttbisABBM3_FD_assembly",
#   "AttbisABBM3_FD_reads",
#   "AttcapACBM1_FD_assembly",
#   "AttcapACBM1_FD_reads",
#   "AttcapACBM2_FD_assembly",
#   "AttcapACBM2_FD_reads",
#   "AttcapACBM3_FD_assembly",
#   "AttcapACBM3_FD_reads",
# "AttlaeALBM1_FD_assembly",
# "AttlaeALBM1_FD_reads",
# "AttlaeALBM2_FD_assembly",
# "AttlaeALBM2_FD_reads",
# "AttlaeALBM3_FD_assembly",
# "AttlaeALBM3_FD_reads",
# "AttsexASBM1_FD_assembly",
# "AttsexASBM1_FD_reads",
# "AttsexASBM2_FD_assembly",
# "AttsexASBM2_FD_reads",
# "AttsexASBM3_FD_assembly",
# "AttsexASBM3_FD_reads",
# "FG1_reads",
# "FG2_reads",
# "FG3_reads")

#gaps_col <- seq(2, 24, by = 2)

pheatmap(filt_data[,cols_ord], 
         cellheight = 10, 
         cellwidth = 10,
         filename = 'output/detection_heatmap.pdf',
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = gaps_col)



pheatmap(filt_data[,], 
         cellheight = 10, 
         cellwidth = 10,
         filename = 'output/detection_heatmap_clustered.pdf',
         cluster_rows = FALSE, annotation_row = db_map)


pheatmap(filt_data[,cols_ord], 
         cellheight = 10, 
         cellwidth = 10,
         filename = 'output/detection_heatmap.png',
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         gaps_col = gaps_col)

pheatmap(filt_data[,], 
         cellheight = 10, 
         cellwidth = 10,
         filename = 'output/detection_heatmap_clustered.png',
         cluster_rows = FALSE)

save(id_map, file = 'id_map.Rdata')











filter(filt_data, row_id == "Leucoagaricus gongylophorus:GCA_022457215.1")

filter(filt_data, row_id == "Leucoagaricus gongylophorus:GCA_022457215.1") %>%
  ggplot(aes(y = intersect_bp, x = id)) + geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))








