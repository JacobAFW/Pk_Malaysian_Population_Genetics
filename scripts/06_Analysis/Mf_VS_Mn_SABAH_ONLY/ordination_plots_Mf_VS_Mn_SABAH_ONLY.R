# Load packages
library(tidyverse)
library(ape)
library(ggtree)

# PCA

## Read in data
pca <- read_table("Pk.eigenvec", col_names=F) %>%
    select(-1) %>%
    mutate(X2 = str_remove(X2, "_DK.*")) %>%
    rename("sampleid" = "X2") %>%
    rename_at(vars(starts_with("X")), ~str_replace(., "X.*", paste0("PC", seq_along(.))))

eigenval <- scan("Pk.eigenval")

## Add metadata 
pca <- pca %>%
    left_join(read_csv("Pk.csv", col_names = FALSE) %>% 
        select(1, 6, 7) %>%
        rename(sampleid = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>%
        mutate(sampleid = str_remove(sampleid, "_DK.*"))) %>%
        mutate(Cluster = ifelse(grepl("PK_SB_DNA_011", sampleid) | # does the ID match any of these
                                grepl("PK_SB_DNA_091", sampleid) | 
                                grepl("PK_SB_DNA_043", sampleid) |
                                grepl("PK_SB_DNA_016", sampleid) |  
                                grepl("PK_SB_DNA_092", sampleid) | 
                                grepl("PK_SB_DNA_030", sampleid) | 
                                grepl("PK_SB_DNA_093", sampleid) | 
                                grepl("PK_SB_DNA_042", sampleid) | 
                                grepl("PK_SB_DNA_063", sampleid) |
                                grepl("PK_SB_DNA_028", sampleid) , "Mn", "Mf"))


## Percentage variance explained by each PC
pve_plot <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100) %>%
    ggplot(aes(PC, pve)) + 
    geom_bar(stat = "identity") + 
    ylab("Percentage variance explained") + 
    theme_light()

ggsave("Sabah_only_Mn_vs_Mf/Percentage_variance_explained.png", dpi=600, pve_plot)


## PCA - Mf vs Mn
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Cluster)) + 
    geom_point(size = 3) +
    scale_color_manual(values = c("#440154FF", "#39568CFF")) +
    coord_equal() + 
    theme_light() + 
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
    xlim(-0.5, 0.5)

ggsave("Sabah_only_Mn_vs_Mf/PCA.png", dpi = 600, pca_plot)


# MDS - Mf vs Mn

## Read in data
mds <- read_table("Pk.mds", col_names=T) %>%
    select(-c("FID", "X14")) %>%
    mutate(IID = str_remove(IID, "_DK.*")) %>%
    rename("sampleid" = "IID") %>%
    rename_at(vars(starts_with("C")), ~str_replace(., "C", "MDS"))

## Add metadata 
mds <- mds %>%
    left_join(read_csv("Pk.csv", col_names = FALSE) %>% 
        select(1, 6, 7) %>%
        rename(sampleid = X1) %>%
        rename(Location = X6) %>%
        rename(Cluster = X7) %>%
        mutate(sampleid = str_remove(sampleid, "_DK.*"))) %>%
        mutate(Cluster = ifelse(grepl("PK_SB_DNA_011", sampleid) | # does the ID match any of these
                                grepl("PK_SB_DNA_091", sampleid) | 
                                grepl("PK_SB_DNA_043", sampleid) |
                                grepl("PK_SB_DNA_016", sampleid) |  
                                grepl("PK_SB_DNA_092", sampleid) | 
                                grepl("PK_SB_DNA_030", sampleid) | 
                                grepl("PK_SB_DNA_093", sampleid) | 
                                grepl("PK_SB_DNA_042", sampleid) | 
                                grepl("PK_SB_DNA_063", sampleid) |
                                grepl("PK_SB_DNA_028", sampleid) , "Mn", "Mf"))

## Plot MDS
mds_plot <- ggplot(mds, aes(MDS1, MDS2, colour = Cluster)) + 
    geom_point(size = 3) +
    scale_color_manual(values = c("#440154FF", "#39568CFF")) +
    coord_equal() + 
    theme_light()

ggsave("Sabah_only_Mn_vs_Mf/MDS.png", dpi=600, mds_plot)

# Neighbour-joining tree

## Create distance matrix from PLINK data and build NJT
NJT_ID <- read_table("Pk.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DKD.*")) 

NJT_matrix <- read_table("Pk.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

## Plot tree
# options(ignore.negative.edge=TRUE)

# Read in initial metadata
NJT_metadata <- read_csv("Pk.csv", col_names = FALSE) %>%
    select(1, 6, 7) %>%
    rename(Sample = X1) %>%
    rename(Location = X6) %>%
    rename(Cluster = X7) %>%
    mutate(Sample = str_remove(Sample, "_DKD.*")) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(Region = ifelse(Location == "Sabah" | Location == "Betong" | Location == "Kapit" | Location == "Sarikei", "Borneo", "Peninsular"))  %>%
    mutate(Cluster = as.character(Cluster)) %>%
    mutate(Cluster = ifelse(grepl("PK_SB_DNA_011", Sample) | # does the ID match any of these
                            grepl("PK_SB_DNA_091", Sample) | 
                            grepl("PK_SB_DNA_043", Sample) |
                            grepl("PK_SB_DNA_016", Sample) |  
                            grepl("PK_SB_DNA_092", Sample) | 
                            grepl("PK_SB_DNA_030", Sample) | 
                            grepl("PK_SB_DNA_093", Sample) | 
                            grepl("PK_SB_DNA_042", Sample) | 
                            grepl("PK_SB_DNA_063", Sample), "Mn", .$Cluster)) %>%
    mutate(Cluster = ifelse(Cluster == "Sabah", "Mf", .$Cluster)) # convert the remaining Sabah samples to Mf

options(ignore.negative.edge=TRUE)

# Cluster - Mn vs Mf
NJT_tree_plot <- ggtree(NJT_tree, layout="circular", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF"))
    
ggsave("Sabah_only_Mn_vs_Mf/NJT_tree_cluster_rooted.png", dpi = 600, height = 8, width = 16, NJT_tree_plot)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_manual(values = c("#440154FF", "#39568CFF"))
    
ggsave("Sabah_only_Mn_vs_Mf/NJT_tree_cluster_unrooted.png", dpi = 600, height = 15, width = 15, limitsize = FALSE, NJT_tree_plot)

#Labels
NJT_tree_plot <- ggtree(NJT_tree, layout="circular", size = 0.5, aes(colour = Cluster)) %<+% NJT_metadata +
   theme(legend.position = "right", 
   legend.title = element_blank(), 
   legend.key = element_blank()) +
   geom_tiplab(aes(label = Region), size = 2) +
   scale_color_manual(values = c("#440154FF", "#39568CFF"))
    
ggsave("Sabah_only_Mn_vs_Mf/NJT_tree_labelled_rooted.png", dpi = 600, height = 20, width = 20, limitsize = FALSE, NJT_tree_plot)