library(tidyverse)

## NEW DATA MUNGE after defense

load('data/rdata/sourmash_taxa_data_all.rdata')

source('scripts/samples_order_2.R') ## read in cols_ord
cols_ord <- rev(cols_ord)
contam <- scan('data/reagent_lab_contaminating_genera.txt', what = 'character',
               comment.char = "#")

read_as_data <- data_in
## from read_as I need to summarize for the genus, the sample_id and the detection 
## from the 'intersect_bp' section as well as the seq_type

read_as_data2 <- read_as_data %>% select(id, seq_type, genus, intersect_bp) %>%
  group_by(id, seq_type, genus) %>%
  mutate(detected = ifelse(sum(intersect_bp > 0), yes = TRUE, no = FALSE)) %>%
  ungroup() %>%
  select(-intersect_bp)

head(read_as_data2)

load('data/marker_data_counts.Rdata') ## read in "marker_data"

## from marker_data I need to summarize the dection from the 'count' column and preserve
## the sample_id column 'id' as well as the seq type and the genus

marker_data2 <- marker_data %>% 
  rename('seq_type' = seq) %>%
  group_by(id, seq_type, genus) %>%
  mutate(detected = ifelse(sum(count >= 10), yes = TRUE, no = FALSE)) %>%
  ungroup() %>%
  select(-count, -acc)


## combine the datasets

genera_int <- read_csv('output/genera_of_interest.csv') %>%
  filter(genus %in% c(marker_data2$genus, read_as_data2$genus)) 

genera_int_join = genera_int %>%   ## only keep what we have
  mutate(detected = FALSE) %>% 
  mutate(id = list(cols_ord)) %>% ## add in all samples
  mutate(seq_type = list(all_seq_types)) %>%
  unnest(id) %>%
  unnest(seq_type) %>% 
  rename('sample_id'= id) %>%
  select(sample_id, seq_type, genus, detected)

annotation = read_csv('tables/figure_2_colors.csv') %>%
  mutate(id = sub("_.*$", "", id)) %>%
  mutate(sample_id = sub("_read", "", sample_id)) %>%
  mutate(sample_id = sub("_assembly", "", sample_id)) %>% 
  distinct()

all_data <- rbind(read_as_data2, marker_data2)  %>%
  left_join(., annotation, by = 'id') %>%
  filter(!genus %in% contam,
         sample_id %in% cols_ord) %>%
  select(sample_id, seq_type, genus, detected) %>%
  bind_rows(., genera_int_join)


all_seq_types <- all_data$seq_type %>% unique()


data_sub <- all_data %>% 
  group_by(genus, sample_id, seq_type) %>%
  mutate(detected = ifelse(any(detected), yes = sum(na.omit(c(detected, detected))), no = 0)) %>%
  mutate(detected = ifelse(detected > 0, yes = 1, no = 0)) %>%
  ungroup() %>%
  mutate(val = 1) %>%
  distinct() %>%
  left_join(., genera_int, by = 'genus') %>%
  unite(col = 'legend', seq_type, detected)

test <- filter(data_sub, sample_id %in% c("Apterostigma_dentigerum",
                                          "Atta_sexdens_3"), 
               genus %in% c('Enterobacter',
                            'Phialophora'))

View(test)





order_init <- data_sub %>% 
  select(genus, sample_id, taxonomy) %>%
  distinct() %>%
  select(-sample_id) %>%
  group_by(taxonomy) %>%
  mutate(n_genera = length(unique(genus))) %>%
  ungroup() %>%
  distinct()

order <- c(
  "Tremellomycetes", 
  "Saccharomycetes", 
  "Dothideomycetes", 
  "Agaricomycetes", 
  "Sordariomycetes", 
  "Mucoromycetes", 
  "Microbotryomycetes", 
  "Eurotiomycetes", 
  "Gammaproteobacteria", 
  "Actinomycetia", 
  "Bacteroidia",
  "Betaproteobacteria", 
  "Bacilli", 
  "Alphaproteobacteria"
)

x <- order[1]

reorder_genus <- lapply(order, function(x){
  order_init %>% filter(taxonomy == x)
}) %>% bind_rows()

## pick up here
fill_cols = c("read_1","assembly_1","16SrRNA_0",
"16SrRNA_1","18SrRNA_1","18SrRNA_0","28SrRNA_1","28SrRNA_0",
"ITS_0",
"ITS_1",
"assembly_0",
"read_0")

data_sub$sample_id  <- factor(data_sub$sample_id, levels = cols_ord)
data_sub$genus      <- factor(data_sub$genus, levels = reorder_genus$genus)

data_sub$seq_type      <- factor(data_sub$seq_type, levels = all_seq_types)



## y = sample_id
## x = genus
## fill = detected/seq_type


ggplot(test, aes(x=detected, y = val, fill = seq_type)) +
  geom_bar(position="fill", stat="identity", width = .1) + 
  facet_grid(sample_id ~ genus) + 
#  coord_polar("y") + 
  theme_classic() + 
  labs(x = 'Sample', y = 'Genus') +
  theme(legend.position = 'bottom',
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) 


ggsave(plot = p, filename = '_figures_main/Figure3/figure3_v2.pdf', 
       width = 45, height = 45,
       limitsize = FALSE)


############ testing ############
mydata <- data.frame(side1=rep(LETTERS[1:3],3,each=9),
                     side2=rep(LETTERS[1:3],9,each=3),
                     widget=rep(c("X","Y","Z"),9*3),
                     val=runif(9*3),
                     strength=rep(c(1,2,3),3,each=3))



ggplot(mydata, aes(x=strength/2, y = val, fill = widget, width = strength)) +
  geom_bar(position="fill", stat="identity") + 
  facet_grid(side1 ~ side2) + 
  coord_polar("y")





