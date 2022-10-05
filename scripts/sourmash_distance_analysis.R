### distances


# The `compare` subcommand compares one or more signatures (created with
#                                                           `sketch`) using estimated Jaccard index [1] or (if signatures are
#                                                                                                           created with `-p abund`) the angular similarity [2]).
# 
# The default output is a text display of a similarity matrix where each
# entry `[i, j]` contains the estimated Jaccard index between input
# signature `i` and input signature `j`.  The output matrix can be saved
# to a file with `--output` and used with the `sourmash plot` subcommand
# (or loaded with `numpy.load(...)`.  Using `--csv` will output a CSV
#   file that can be loaded into other languages than Python, such as R.


library(tidyverse)
library(pheatmap)

#load('id_map.Rdata')

#id_map_l <- id_map$id

#names(id_map_l) <- id_map$query_filename

df <- read_csv('sourmash/compare/reads_and_assemblies_all.csv') 


df$row <- colnames(df)

df_dist <- df %>% column_to_rownames('row')


pheatmap(df_dist, 
         cellwidth = 8,
         cellheight = 8,
         filename = 'output/distance_heatmap.pdf')
pheatmap(df_dist, 
         cellwidth = 8,
         cellheight = 8,
         filename = 'output/distance_heatmap.png')
