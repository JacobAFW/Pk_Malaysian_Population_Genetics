# Load packages
library(tidyverse)
library(ape)
library(ggtree)

# Neighbour-joining tree

## Create distance matrix from PLINK data and build NJT
NJT_ID <- read_table("DR_01.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("DR_01.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

# Read in initial metadata
NJT_metadata <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/outputs/05_Analyses/Pk.csv", col_names = FALSE) %>% 
        select(1, 6, 7) %>%
        rename(Sample = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>% 
        as.data.frame() 

NJT_metadata <- NJT_metadata %>% 
    mutate(Cluster = ifelse(grepl("PK_SB_DNA_011", Sample) | # does the ID match any of these
                                grepl("PK_SB_DNA_091", Sample) | 
                                grepl("PK_SB_DNA_043", Sample) |
                                grepl("PK_SB_DNA_016", Sample) |  
                                grepl("PK_SB_DNA_092", Sample) | 
                                grepl("PK_SB_DNA_030", Sample) | 
                                grepl("PK_SB_DNA_093", Sample) | 
                                grepl("PK_SB_DNA_042", Sample) | 
                                grepl("PK_SB_DNA_063", Sample), "Mn", .$Cluster)) %>% # if not, just use values from X7 - clusters and Sabah
    mutate(Cluster = ifelse(Cluster == "Sabah", "Mf", .$Cluster)) %>% # if its the remaining Sabah samples, make them Mn, else keep them the same
    mutate_all(~str_remove(., "_DKD.*")) 

# Plot tree
## options(ignore.negative.edge=TRUE)

options(ignore.negative.edge=TRUE)

# 01
# Cluster 
NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF"))
    
ggsave("NJT_unrooted_01.png", dpi = 600, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

#Labels
NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
    legend.title = element_blank(), 
    legend.key = element_blank()) +
    geom_tiplab(size = 2) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF"))

ggsave("NJT_labelled_unrooted_01.png", dpi = 600, height = 20, width = 20, limitsize = FALSE, NJT_tree_plot)

# 05
NJT_ID <- read_table("DR_05.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("DR_05.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)
NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF"))
    
ggsave("NJT_unrooted_05.png", dpi = 600, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

# 10
NJT_ID <- read_table("DR_10.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("DR_10.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF"))
    
ggsave("NJT_unrooted_10.png", dpi = 600, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

# 13
NJT_ID <- read_table("DR_13.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("DR_13.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF"))
    
ggsave("NJT_unrooted_13.png", dpi = 600, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

# 14.1
NJT_ID <- read_table("DR_14.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("DR_14.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF"))
    
ggsave("NJT_unrooted_14.1.png", dpi = 600, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

# 14.2
NJT_ID <- read_table("DR_14_2.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("DR_14_2.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF"))
    
ggsave("NJT_unrooted_14.2.png", dpi = 600, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)
